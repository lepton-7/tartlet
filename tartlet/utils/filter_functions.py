from tart.utils.activity_inference import Candidate


def default_check(cand: Candidate):
    issues = []
    warns = []

    if abs(cand.from_switch_end_relative) > 0.3:
        issues.append(f"Too far from end (abs({cand.from_switch_end_relative})>0.3)")
    if cand.rel_cov_delta > -0.2:
        issues.append(f"Small drop ({cand.rel_cov_delta}>-0.2)")
    if cand.symks_pval > 0.05:
        warns.append("Symmetry check failed (>0.05)")

    if len(issues) > 0:
        cand.note = ";".join(issues)
        return "fail"

    else:
        cand.note = ";".join(warns)

    return "pass"
