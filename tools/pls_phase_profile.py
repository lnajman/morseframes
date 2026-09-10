"""Validate optional fine PLS timings without double-counting parent phases."""
import math

GROUPS = {
    "setup": ("output_init", "vertex_order", "executor_init", "storage_init", "owner_keys", "partition"),
    "local": ("schedule", "execution"),
    "cleanup": ("events_index_cleanup", "keys_cleanup", "membership_cleanup", "executor_cleanup", "vertices_cleanup"),
}
FIELDS = {name for group in GROUPS.values() for name in group} | {"cleanup"}


def validate(profile, setup, local, replay, gradient):
    if set(profile) != FIELDS or any(not math.isfinite(v) or v < 0 for v in profile.values()):
        raise ValueError("Invalid fine PLS profile fields")
    result = dict(profile)
    for group, parent in [("setup", setup), ("local", local), ("cleanup", profile["cleanup"])]:
        residual = parent - sum(profile[k] for k in GROUPS[group])
        if residual < -max(1e-12, parent * 1e-9):
            raise ValueError("Fine PLS phases exceed their parent")
        if group == "cleanup" and abs(residual) > max(1e-12, parent * 1e-9):
            raise ValueError("Cleanup checkpoints do not partition cleanup")
        result[group + "_unattributed"] = max(0., residual)
    residual = gradient - setup - local - replay - profile["cleanup"]
    if residual < -max(1e-12, gradient * 1e-9):
        raise ValueError("PLS cleanup exceeds gradient remainder")
    result["gradient_unattributed"] = max(0., residual)
    return result
