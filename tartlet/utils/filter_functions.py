from tart.utils.activity_inference import Candidate


class DefaultThresholds:
    relative_size_bounds = 0.3
    relative_cov_delta = -0.2
    symmetric_peaks_pval = 0.05

    def __init__(self) -> None:
        pass
    
    @staticmethod
    def check(cand: Candidate):
        issues = []
        warns = []

        if abs(cand.from_switch_end_relative) > DefaultThresholds.relative_size_bounds:
            issues.append(
                f"Too far from end (abs({cand.from_switch_end_relative})>{DefaultThresholds.relative_size_bounds})"
            )
        if cand.rel_cov_delta > DefaultThresholds.relative_cov_delta:
            issues.append(
                f"Small drop ({cand.rel_cov_delta}>{DefaultThresholds.relative_cov_delta})"
            )
        if cand.symks_pval > DefaultThresholds.symmetric_peaks_pval:
            warns.append(
                f"Symmetry check failed (>{DefaultThresholds.symmetric_peaks_pval})"
            )

        if len(issues) > 0:
            cand.note = ";".join(issues)
            return "fail"

        else:
            cand.note = ";".join(warns)

        return "pass"
