"""Background worker thread that runs an algorithm and reports progress."""

import threading

from gui.registries import make_algorithm


class Worker(threading.Thread):
    def __init__(self, msg_q, problem, algo_key, pop_size, max_nfe, update_freq):
        super().__init__(daemon=True)
        self.q = msg_q
        self.problem = problem
        self.algo_key = algo_key
        self.pop_size = pop_size
        self.max_nfe = max_nfe
        self.update_freq = update_freq
        self._stop_event = threading.Event()

    def request_stop(self):
        self._stop_event.set()

    def run(self):
        try:
            algo = make_algorithm(self.algo_key, self.problem, self.pop_size)
            state = {"last": 0}

            def cb(a):
                if self._stop_event.is_set():
                    raise InterruptedError
                if a.nfe - state["last"] >= self.update_freq:
                    state["last"] = a.nfe
                    snap = [list(s.objectives[:]) for s in a.result]
                    self.q.put({"t": "update", "nfe": a.nfe, "result": snap})

            algo.run(self.max_nfe, callback=cb)
            final = [list(s.objectives[:]) for s in algo.result]
            self.q.put({"t": "done", "nfe": algo.nfe, "result": final})

        except InterruptedError:
            self.q.put({"t": "stopped"})
        except Exception as e:
            import traceback
            self.q.put({"t": "error", "msg": str(e), "tb": traceback.format_exc()})
