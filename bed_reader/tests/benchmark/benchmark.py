if __name__ == "__main__":

    import numpy as np
    from pathlib import Path
    from bed_reader import to_bed
    import time
    import pandas as pd
    import matplotlib.pyplot as plt

    if True:
        ssd_path = Path("m:/deldir/bench")
        hdd_path = Path("e:/deldir/bench")
    else:
        ssd_path = Path("/mnt/m/deldir/bench")
        hdd_path = Path("/mnt/e/deldir/bench")

    def test_writes(
        iid_count, sid_count, num_threads, drive, include_error, version_list
    ):
        if drive == "ssd":
            path = ssd_path
        elif drive == "hdd":
            path = hdd_path
        else:
            raise ValueError(f"drive must be 'ssd' or 'hdd', not '{drive}'")

        output_file = path / f"{iid_count}x{sid_count}.bed"

        val = np.full((iid_count, sid_count), 1, dtype=np.int8)
        if include_error:
            val[iid_count // 2, sid_count // 2] = 22

        val_size = float(iid_count) * sid_count
        if val_size > 9_200_000_000:
            raise ValueError(f"val_size {val_size} is too large")

        result_list = []
        for version in version_list:
            start = time.time()
            to_bed(output_file, val, num_threads=num_threads, version=version)
            delta = time.time() - start
            result = [
                iid_count,
                sid_count,
                num_threads,
                drive,
                include_error,
                version,
                True,
                "i8",
                "benchmark.py",
                val_size,
                round(delta, 4),
                round(val_size / delta),
                round(val_size / delta / num_threads),
            ]
            print(result)
            result_list.append(result)

        result_df = pd.DataFrame(
            result_list,
            columns=[
                "iid_count",
                "sid_count",
                "num_threads",
                "drive",
                "include_error",
                "version",
                "release",
                "dtype",
                "source",
                "val size",
                "time",
                "val per second",
                "per thread",
            ],
        )
        return result_df

    def meta_test(
        iid_count,
        sid_start=5,
        sid_end=None,
        point_count=30,
        drive_list=["ssd"],
        plot_index=0,
    ):
        # 5K vs 50K
        # 50K vs 50K
        # 500K vs 5K
        if sid_end is None:
            if iid_count <= 50_000:
                sid_end = 50_000
            else:
                sid_end = 5_000

        result = []
        for sid_count in np.logspace(
            np.log10(sid_start), np.log10(sid_end), point_count, base=10, dtype=int
        ):
            for drive in drive_list:
                for num_threads in [1, 12]:
                    result.append(
                        test_writes(
                            iid_count, sid_count, num_threads, drive, False, [3, 4]
                        )
                    )
        df = pd.concat(result)
        df2 = df.pivot(
            index="sid_count",
            columns=["iid_count", "drive", "num_threads", "version"],
            values="val per second",
        )
        df2.plot(marker=".", logx=True)
        plt.savefig(
            ssd_path
            / "plots"
            / f"plot{plot_index},{iid_count},{'_'.join(drive_list)}.png"
        )
        # plt.show()
        return df


plot_count = 0
for drive in ["ssd", "hdd"]:
    for iid_count in [5_000, 50_000, 500_000]:
        meta_test(
            iid_count,
            drive_list=[drive],
            plot_index=plot_count,
            sid_end=500,
            point_count=10,
        )
        plot_count += 1
