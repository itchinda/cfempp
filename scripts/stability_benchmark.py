import subprocess
import re
import statistics
import argparse
import matplotlib.pyplot as plt


# ---------------------------
# CLI
# ---------------------------
parser = argparse.ArgumentParser(description="CFEM benchmark runner")
parser.add_argument("--runs", type=int, default=10, help="Number of runs")
parser.add_argument("--exe", type=str, default="./apps/cfem_app_opt", help="CFEM executable")
args = parser.parse_args()

N_RUNS = args.runs
EXECUTABLE = args.exe


# ---------------------------
# PARSING FUNCTION
# ---------------------------
def extract_real_time(output: str) -> float:
    """
    Supports both:
    - GNU time format (-p): real 56.61
    - formatted format: real 0m56.608s
    """

    # Case 1: -p format (most common)
    match = re.search(r"real\s+([\d\.]+)", output)
    if match:
        return float(match.group(1))

    # Case 2: classic formatted output
    match = re.search(r"real\s+(\d+)m([\d\.]+)s", output)
    if match:
        minutes = float(match.group(1))
        seconds = float(match.group(2))
        return 60 * minutes + seconds

    # Debug aid
    raise ValueError("Could not parse time output:\n" + output)
# ---------------------------
# RUN BENCHMARK
# ---------------------------
times = []

for i in range(N_RUNS):
    print(f"Run {i+1}/{N_RUNS}")

    result = subprocess.run(
        f"/usr/bin/time -p {EXECUTABLE}",
        shell=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True
    )

    output = result.stderr

    t = extract_real_time(output)
    times.append(t)

    print(f"  -> {t:.3f} s")


# ---------------------------
# STATISTICS
# ---------------------------
avg = statistics.mean(times)
std = statistics.stdev(times) if len(times) > 1 else 0.0
tmin = min(times)
tmax = max(times)


# ---------------------------
# PLOT
# ---------------------------
plt.figure(figsize=(10, 5))

plt.plot(times, marker='o', label="Run time")
plt.axhline(avg, linestyle='--', label=f"Average = {avg:.2f}s")

plt.title("CFEM Execution Time Stability")
plt.xlabel("Run index")
plt.ylabel("Time (s)")
plt.legend()

text = (
    f"avg = {avg:.3f} s\n"
    f"std = {std:.3f} s\n"
    f"min = {tmin:.3f} s\n"
    f"max = {tmax:.3f} s"
)

plt.gca().text(
    0.02, 0.95, text,
    transform=plt.gca().transAxes,
    verticalalignment='top',
    bbox=dict(boxstyle="round", facecolor="white", alpha=0.8)
)

plt.grid(True)
plt.tight_layout()
plt.savefig('stability_benchmark.png', dpi=300)
plt.show()
