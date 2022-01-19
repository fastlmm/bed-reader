import numpy as np
from pathlib import Path
from bed_reader import to_bed
import time
import pandas as pd
import matplotlib.pyplot as plt

ssd_path = Path(f"m:/deldir/bench")
hdd_path = Path(f"e:/deldir/bench")


def test_writes(iid_count, sid_count, num_threads, drive, include_error, algo):
    if drive == "ssd":
        path = ssd_path
    elif drive == "hdd":
        path = hdd_path
    else:
        raise ValueError(f"drive must be 'ssd' or 'hdd', not '{drive}'")

    output_file = path / f"{iid_count}x{sid_count}.bed"

    val = np.full((iid_count, sid_count), 1, dtype=np.float64)
    if include_error:
        val[iid_count // 2, sid_count // 2] = 22

    val_size = iid_count * sid_count
    if val_size > 10_000_000_000:
        raise ValueError(f"val_size {val_size} is too large")

    result_list = []
    for version in algo:
        start = time.time()
        to_bed(output_file, val, version=version, num_threads=num_threads)
        delta = time.time() - start
        num_threads_adjusted = num_threads if version > 0 else 1
        result = [
            iid_count,
            sid_count,
            num_threads,
            drive,
            include_error,
            version,
            True,
            "benchmark.py",
            val_size,
            round(delta, 4),
            round(val_size / delta),
            round(val_size / delta / num_threads_adjusted),
        ]
        result_list.append(result)
        print(result)

    result_df = pd.DataFrame(
        result_list,
        columns=[
            "iid_count",
            "sid_count",
            "num_threads",
            "drive",
            "include_error",
            "algorithm",
            "release",
            "source",
            "val size",
            "time",
            "val per second",
            "per thread",
        ],
    )
    return result_df


if False:
    test_writes(1_000, 1_000, 1, "ssd", False, [0, 1, 1, 0])
    test_writes(1_000, 1_000, 1, "hdd", False, [0, 1, 1, 0])

if False:
    result = []
    for sid_count in np.logspace(np.log10(5), np.log10(62_000), 25, base=10, dtype=int):
        iid_count = 10_000
        for drive in ["ssd"]:
            for num_threads in [12]:
                result.append(
                    test_writes(iid_count, sid_count, num_threads, drive, False, [0, 1])
                )
    df = pd.concat(result)
    df2 = df.pivot(
        index="sid_count",
        columns=["algorithm", "drive", "num_threads"],
        values="val per second",
    )
    df2.plot(marker=".", logx=True)
    plt.show()

if True:
    result = []
    for sid_count in np.logspace(np.log10(5), np.log10(35_000), 25, base=10, dtype=int):
        iid_count = 50_000
        for drive in ["ssd"]:
            for num_threads in [12]:
                result.append(
                    test_writes(iid_count, sid_count, num_threads, drive, False, [1])
                )
    df = pd.concat(result)
    df2 = df.pivot(
        index="sid_count",
        columns=["algorithm", "drive", "num_threads"],
        values="val per second",
    )
    df2.plot(marker=".", logx=True)
    plt.show()
