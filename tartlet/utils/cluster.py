import numpy as np
import pandas as pd

from random import sample
from tart.utils.utils import print
from scipy.stats import mannwhitneyu, levene
from scipy.cluster.hierarchy import fclusterdata


class Cluster:
    def __init__(
        self,
        peak_log: pd.DataFrame,
        locus_column: str = "rowid",
        feature_dimension_columns: list[str] = ["from_riboswith_end_relative"],
    ) -> None:

        self.rowid = locus_column
        self.feature_dims = feature_dimension_columns

        self.statdim = "coverage_delta_stable_relative"
        self.posdim = "from_riboswith_end_relative"

        self.df = self.__cluster(peak_log)
        self.stats_df = self.__compute_stats()

    def get(self):
        return (self.df, self.stats_df)

    def __cluster(self, df) -> pd.DataFrame:
        tdfs: list[pd.DataFrame] = []

        try:
            _ = df[self.rowid]
        except KeyError:
            print(f"Df empty for clustering. Returning unmodified argument.")
            return df

        for rowid in pd.unique(df[self.rowid]):
            try:
                tdf = df[df[self.rowid] == rowid].dropna(subset=self.feature_dims)
            except:  # FIXME: Fix this properly
                print(f"Skipping id: {rowid}.")
                continue

            # FIXME: Make this somehow compatible with arbitrary dims and col names
            relpos = np.array(tdf[self.posdim])

            # Need to reshape the data for clustering
            relpos_reshaped = np.reshape(relpos, (len(relpos), len(self.feature_dims)))

            # This is the fully parameterised clustering function that returns labels.
            # A pretty round-way of abstracting away clustering parameters that will probably
            # slow everything down ¯\_(ツ)_/¯. Hopefully this will allow easy user customisability
            # in the future
            clfunc = self.__default_cluster(relpos_reshaped)

            try:
                tdf["cluster"] = clfunc()
            except:  # FIXME: Figure out what makes clustering break
                print(f"Clustering broke. Skipping.")
                continue

            tdfs.append(tdf)

        redf = pd.concat(tdfs, ignore_index=True)
        return df.merge(redf, how="left")

    def __compute_stats(self) -> pd.DataFrame:
        skel_list = []

        try:
            _ = self.df[self.rowid]
        except KeyError:
            print(f"Df empty for stats. Returning empty")
            return self.df

        # Make the skeleton dataframe containing rowid and cluster columns
        for rowid in pd.unique(self.df[self.rowid]):
            tdf = self.df[self.df[self.rowid] == rowid].dropna(subset=self.feature_dims)
            for cl in pd.unique(tdf["cluster"]):
                skel_list.append({self.rowid: rowid, "cluster": cl})

        df = pd.DataFrame(skel_list)

        pos_mean_arr = []
        pos_var_arr = []
        delta_mean_arr = []
        delta_var_arr = []
        delta_meanp_arr = []
        delta_varp_arr = []
        exset_delta_mean_arr = []
        exset_delta_variance_arr = []

        # Iterate over each rowid and cluster combination
        for i, row in df.iterrows():
            # compute mean and variance pvals for each combination
            rowid = row[self.rowid]
            cl = row["cluster"]

            # Set of peaks in cluster
            peakset = self.df[
                (self.df[self.rowid] == rowid) & (self.df["cluster"] == cl)
            ].dropna(subset=self.feature_dims)

            # Set of cluster exclusive peaks
            exset = self.df[
                (self.df[self.rowid] == rowid) & (self.df["cluster"] != cl)
            ].dropna(subset=self.feature_dims)

            # Test if peakset mean delta is < exset mean delta
            try:
                meanp = mannwhitneyu(
                    x=peakset[self.statdim], y=exset[self.statdim], alternative="less"
                ).pvalue
            except ValueError:
                meanp = 1.0

            # Test against null -> peakset delta variance is ~= exset delta variance
            try:
                varp = self.__var_test(peakset[self.statdim], exset[self.statdim])
            except ValueError:
                varp = 1.0

            pos_mean_arr.append(np.mean(peakset[self.posdim]))
            pos_var_arr.append(np.var(peakset[self.posdim]))
            delta_mean_arr.append(np.mean(peakset[self.statdim]))
            delta_var_arr.append(np.var(peakset[self.statdim]))
            delta_meanp_arr.append(meanp)
            delta_varp_arr.append(varp)
            exset_delta_mean_arr.append(np.mean(exset[self.statdim]))
            exset_delta_variance_arr.append(np.var(exset[self.statdim]))

        df["pos_mean"] = pos_mean_arr
        df["pos_variance"] = pos_var_arr
        df["delta_mean"] = delta_mean_arr
        df["noiseset_delta_mean"] = exset_delta_mean_arr
        df["delta_mean_pval"] = delta_meanp_arr
        df["delta_variance"] = delta_var_arr
        df["noiseset_delta_variance"] = exset_delta_variance_arr
        df["delta_variance_pval"] = delta_varp_arr

        return df

    def __var_test(self, cut: pd.Series, exset: pd.Series):
        n = len(cut)
        exsize = len(exset)

        # For the resampling to avoid a div by 0 error
        adj = 1 if n > exsize else 0

        # The CUT isn't large enough to begin with
        if n < 3:
            return 1.0

        exshuff = list(exset)

        # shuffle then subset to emulate in-place sampling-without-replacement
        for _ in range(5):
            exshuff = sample(exshuff, exsize)

        # Collect at least rep comparisons to get a distribution of pvals
        rep = 100
        ex_collection = []
        for _ in range((rep // (exsize // n + adj)) + 1):
            ex_collection.extend(
                [exshuff[i : i + n] for i in range(0, exsize, n) if i + n <= exsize]
            )

        # Levene's test on CUT against every sampled exset subset
        plist = np.array(
            map(lambda ex: levene(cut, ex, center="median").pvalue, ex_collection)
        )

        return np.percentile(plist, 50)

    def __default_cluster(self, X: np.ndarray):
        def wrapper():
            return fclusterdata(X, t=0.04, criterion="distance", method="complete")

        return wrapper
