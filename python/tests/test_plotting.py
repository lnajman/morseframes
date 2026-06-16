import importlib.util
import unittest

import morseframes as mp


@unittest.skipUnless(importlib.util.find_spec("plotly") is not None, "Plotly is not importable")
class PlottingTest(unittest.TestCase):
    def test_two_panel_surface_demo_uses_one_sequence(self):
        from morseframes import plotting

        field = plotting.build_noisy_sine_square(size=4, noise=0.0, seed=1)
        self.assertIsNotNone(field.original_values)
        frame = mp.compute_morse_sequence_and_reference_map(field.complex, algorithm="f-max")
        coreferences = mp.compute_coreference_map(
            field.complex,
            frame.sequence,
            algorithm="f-max",
        )

        critical_simplexes = frame.sequence.critical_simplices_as_simplices(field.complex)
        cocycle_supports = plotting.reference_cocycle_edges(
            field.complex,
            frame.sequence,
            frame.references,
        )
        cycle_supports = plotting.coreference_cycle_edges(
            field.complex,
            frame.sequence,
            coreferences,
        )

        self.assertGreater(len(critical_simplexes), 0)
        for critical_edge, support_edges in cocycle_supports.items():
            self.assertEqual(len(critical_edge), 2)
            self.assertIn(critical_edge, critical_simplexes)
            for edge in support_edges:
                self.assertEqual(len(edge), 2)
                self.assertIn(edge, field.complex)
        for critical_edge, support_edges in cycle_supports.items():
            self.assertEqual(len(critical_edge), 2)
            self.assertIn(critical_edge, critical_simplexes)
            for edge in support_edges:
                self.assertEqual(len(edge), 2)
                self.assertIn(edge, field.complex)

        fig = plotting.plot_morse_surface_with_persistence(
            field,
            algorithm="f-max",
            persistence_steps=4,
        )
        payload = fig.to_plotly_json()
        self.assertIn("scene", payload["layout"])
        self.assertEqual(len(payload["layout"]["sliders"]), 2)
        self.assertEqual(payload["layout"]["legend"]["orientation"], "v")
        self.assertEqual(len(payload["layout"]["updatemenus"]), 1)
        self.assertEqual({frame["name"] for frame in payload["frames"]}, {"noisy", "original"})
        self.assertTrue(any(trace["type"] == "mesh3d" for trace in payload["data"]))
        self.assertTrue(any(trace["type"] == "scatter3d" for trace in payload["data"]))
        self.assertTrue(any(trace["type"] == "scatter" for trace in payload["data"]))
        trace_names = {trace.get("name") for trace in payload["data"]}
        self.assertIn("critical 0", trace_names)
        self.assertIn("critical 1", trace_names)
        self.assertIn("critical 2", trace_names)
        self.assertIn("1-cycles", trace_names)
        self.assertIn("1-cocycles", trace_names)


if __name__ == "__main__":
    unittest.main()
