"""A step of a pipeline for the tests of its runner, run as a script.

It takes the command and the package like a command of pkpdsim, then burns CPU
time, starts a process which sleeps, sleeps itself and exits with the code it
is given. It writes its process, the CPUs it sees and the start into a record,
and the end once it ended.
"""

import argparse
import json
import os
import signal
import subprocess
import sys
import time
from pathlib import Path


def main() -> None:
    """Run the step."""
    p = argparse.ArgumentParser()
    p.add_argument("command")
    p.add_argument("package")
    p.add_argument("--record", type=Path, required=True)
    p.add_argument("--cores", type=int, default=None)
    p.add_argument("--seconds", type=float, default=0.0)
    p.add_argument("--burn", type=float, default=0.0)
    p.add_argument("--exit", type=int, default=0)
    p.add_argument("--child", action="store_true")
    p.add_argument("--ignore-interrupt", action="store_true")
    args = p.parse_args()
    if args.ignore_interrupt:
        signal.signal(signal.SIGINT, signal.SIG_IGN)
    record = {
        "pid": os.getpid(),
        "cpus": os.process_cpu_count(),
        "cores": args.cores,
        "start": time.time(),
    }
    if args.child:
        child = subprocess.Popen([sys.executable, "-c", "import time; time.sleep(60)"])
        record["child"] = child.pid
    args.record.write_text(json.dumps(record), encoding="utf-8")
    burn = time.process_time() + args.burn
    while time.process_time() < burn:
        pass
    time.sleep(args.seconds)
    record["end"] = time.time()
    args.record.write_text(json.dumps(record), encoding="utf-8")
    print(f"{args.command} done", flush=True)
    sys.exit(args.exit)


if __name__ == "__main__":
    main()
