"""Exercise the live HTTP boundary, including obsolete slider requests."""

import json
import threading
import time
import unittest
from urllib.error import HTTPError
from urllib.request import Request, urlopen

from playground_server import PlaygroundServer


class PlaygroundHTTPTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.server = PlaygroundServer(("127.0.0.1", 0))
        cls.thread = threading.Thread(target=cls.server.serve_forever, daemon=True)
        cls.thread.start()
        cls.url = f"http://127.0.0.1:{cls.server.server_port}"

    @classmethod
    def tearDownClass(cls):
        cls.server.shutdown()
        cls.server.executor.shutdown(wait=True, cancel_futures=True)
        cls.server.server_close()
        cls.thread.join(timeout=5)

    def request(self, path, payload=None, method=None, headers=None):
        request = Request(
            self.url + path,
            data=None if payload is None else json.dumps(payload).encode(),
            headers=headers or {"Content-Type": "application/json"},
            method=method,
        )
        try:
            with urlopen(request, timeout=10) as response:
                return response.status, json.load(response)
        except HTTPError as error:
            return error.code, json.load(error)

    def job_request(self, client, request_id, **overrides):
        return {
            "diagram": "higgs-sunset",
            "stage": "impred",
            "parameters": {},
            "client_id": client,
            "request_id": request_id,
            **overrides,
        }

    def complete(self, identifier):
        deadline = time.monotonic() + 30
        while time.monotonic() < deadline:
            status, result = self.request("/api/jobs/" + identifier)
            self.assertEqual(status, 200)
            if result["status"] in ("complete", "error", "cancelled"):
                return result
            time.sleep(0.05)
        self.fail("Layout job did not complete")

    def test_configuration_and_real_run(self):
        status, configuration = self.request("/api/config")
        self.assertEqual(status, 200)
        self.assertEqual(configuration, self.server.runner.config())
        self.assertEqual(len(configuration["parameters"]), 13)
        status, response = self.request("/api/run", self.job_request("real-run", 1))
        self.assertEqual(status, 202)
        result = self.complete(response["job_id"])
        self.assertEqual(result["status"], "complete", result)
        stages = result["result"]["case"]["stages"]
        selected = next(s for s in stages if s["id"] == "impred")
        self.assertTrue(selected["available"])
        self.assertTrue(selected["constraints"]["valid"])
        self.assertTrue(selected["constraints"]["straight_externals"])

    def test_late_request_cannot_cancel_newer_job(self):
        status, newer = self.request(
            "/api/run", self.job_request("out-of-order", 10, stage="seed")
        )
        self.assertEqual(status, 202)
        status, older = self.request("/api/run", self.job_request("out-of-order", 9))
        self.assertEqual(status, 400)
        self.assertIn("Superseded", older["error"])
        self.assertEqual(self.complete(newer["job_id"])["status"], "complete")

    def test_cancel_running_job(self):
        # A fresh ImPrEd parameter bypasses cached default geometry.
        status, response = self.request(
            "/api/run",
            self.job_request(
                "cancel", 1, parameters={"impred.target": 3.1}, stage="impred"
            ),
        )
        self.assertEqual(status, 202)
        status, cancelled = self.request(
            "/api/jobs/" + response["job_id"], method="DELETE"
        )
        self.assertEqual(status, 200)
        self.assertEqual(cancelled["status"], "cancelled")
        status, response = self.request(
            "/api/run", self.job_request("cancel", 2, stage="seed")
        )
        self.assertEqual(status, 202)
        self.assertEqual(self.complete(response["job_id"])["status"], "complete")

    def test_invalid_inputs_and_cross_origin(self):
        for parameters in (
            {"unknown": 1},
            {"springs.route-rest-scale": 1},
            {"impred.target": float("nan")},
        ):
            status, _ = self.request(
                "/api/run", self.job_request("invalid", 1, parameters=parameters)
            )
            self.assertEqual(status, 400)
        for stage in ("springs", "anneal"):
            status, _ = self.request(
                "/api/run", self.job_request("removed-stage", 1, stage=stage)
            )
            self.assertEqual(status, 400)
        status, _ = self.request(
            "/api/run",
            self.job_request("foreign", 1),
            headers={
                "Content-Type": "application/json",
                "Origin": "https://example.invalid",
            },
        )
        self.assertEqual(status, 403)
        self.assertEqual(self.request("/api/jobs/missing")[0], 404)


if __name__ == "__main__":
    unittest.main()
