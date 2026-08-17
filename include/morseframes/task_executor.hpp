#pragma once

#include <algorithm>
#include <chrono>
#include <condition_variable>
#include <cstddef>
#include <deque>
#include <functional>
#include <future>
#include <memory>
#include <mutex>
#include <stdexcept>
#include <thread>
#include <type_traits>
#include <utility>
#include <vector>

namespace morseframes {

// Fixed-size executor whose worker budget includes the calling thread. Waiting
// threads cooperatively execute queued work, so a level task can safely wait for
// nested facet tasks without exhausting the pool.
class BoundedTaskExecutor {
 public:
  explicit BoundedTaskExecutor(std::size_t max_workers = 0)
      : worker_count_(resolve_worker_count(max_workers)) {
    workers_.reserve(worker_count_ - 1);
    for (std::size_t index = 1; index < worker_count_; ++index) {
      workers_.emplace_back([this]() { worker_loop(); });
    }
  }

  ~BoundedTaskExecutor() {
    {
      std::lock_guard<std::mutex> lock(mutex_);
      stopping_ = true;
    }
    activity_.notify_all();
    for (auto& worker : workers_) {
      worker.join();
    }
  }

  BoundedTaskExecutor(const BoundedTaskExecutor&) = delete;
  BoundedTaskExecutor& operator=(const BoundedTaskExecutor&) = delete;

  std::size_t worker_count() const { return worker_count_; }

  template <typename Function>
  auto submit(Function&& function)
      -> std::future<std::invoke_result_t<std::decay_t<Function>>> {
    using Result = std::invoke_result_t<std::decay_t<Function>>;
    auto task = std::make_shared<std::packaged_task<Result()>>(
        std::forward<Function>(function));
    auto future = task->get_future();
    {
      std::lock_guard<std::mutex> lock(mutex_);
      if (stopping_) {
        throw std::logic_error("Cannot submit work to a stopped executor.");
      }
      tasks_.emplace_back([this, task]() {
        (*task)();
        activity_.notify_all();
      });
    }
    activity_.notify_one();
    return future;
  }

  template <typename Result>
  Result get(std::future<Result>& future) {
    while (future.wait_for(std::chrono::seconds(0)) !=
           std::future_status::ready) {
      if (try_run_one()) {
        continue;
      }
      std::unique_lock<std::mutex> lock(mutex_);
      activity_.wait(lock, [this, &future]() {
        return stopping_ || !tasks_.empty() ||
               future.wait_for(std::chrono::seconds(0)) ==
                   std::future_status::ready;
      });
    }
    return future.get();
  }

 private:
  static std::size_t resolve_worker_count(std::size_t max_workers) {
    if (max_workers != 0) {
      return std::max<std::size_t>(1, max_workers);
    }
    return std::max<std::size_t>(1, std::thread::hardware_concurrency());
  }

  bool try_run_one() {
    std::function<void()> task;
    {
      std::lock_guard<std::mutex> lock(mutex_);
      if (tasks_.empty()) {
        return false;
      }
      task = std::move(tasks_.front());
      tasks_.pop_front();
    }
    task();
    return true;
  }

  void worker_loop() {
    while (true) {
      std::function<void()> task;
      {
        std::unique_lock<std::mutex> lock(mutex_);
        activity_.wait(lock,
                       [this]() { return stopping_ || !tasks_.empty(); });
        if (stopping_ && tasks_.empty()) {
          return;
        }
        task = std::move(tasks_.front());
        tasks_.pop_front();
      }
      task();
    }
  }

  const std::size_t worker_count_;
  std::vector<std::thread> workers_;
  std::deque<std::function<void()>> tasks_;
  std::mutex mutex_;
  std::condition_variable activity_;
  bool stopping_ = false;
};

}  // namespace morseframes
