from tartlet.utils.activity_inference import Candidate


class DefaultThresholds:
    relative_size_bounds = 0.3
    relative_cov_delta = -0.2
    symmetric_peaks_pval = 0.05
    peak_sig_thresh = 0.05

    def __init__(self) -> None:
        pass

    @staticmethod
    def check(
        cand: Candidate,
        rel_cov_change_sig_thresh: float,
        relative_size_bound_thresh: float 
    ):
        d = DefaultThresholds
        issues = []
        warns = []

        if cand.coverage_drop_pvalue >= d.peak_sig_thresh:
            issues.append(
                f"Drop not significant ({cand.coverage_drop_pvalue:.2f} >= {d.peak_sig_thresh})"
            )

        if abs(cand.from_switch_end_relative) > relative_size_bound_thresh:
            issues.append(
                f"Too far from end (abs({cand.from_switch_end_relative:.2f}) > {relative_size_bound_thresh})"
            )
        if cand.stable_rel_cov_delta > float(rel_cov_change_sig_thresh):
            issues.append(
                f"Small drop ({cand.stable_rel_cov_delta:.2f} > {rel_cov_change_sig_thresh})"
            )

        # if cand.symks_pval > d.symmetric_peaks_pval:
        #     warns.append(
        #         f"Symmetry check failed (>{d.symmetric_peaks_pval})"
        #     )

        if len(issues) > 0:
            cand.note = ";".join(issues)
            return "fail"

        else:
            cand.note = ";".join(warns)

        return "pass"
