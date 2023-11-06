from tart.utils.activity_inference import Candidate


def default_check(cand: Candidate):
    issues = []
    warns = []

    if abs(cand.from_switch_end_relative) > 0.2:
        issues.append("Too far from end")
    if cand.rel_cov_delta > -0.2:
        issues.append("Small drop")
    if cand.symks_pval > 0.05:
        warns.append("Symmetry check failed")

    if len(issues) > 0:
        cand.note = ";".join(issues)
        return "fail"

    else:
        cand.note = ";".join(warns)

    return "pass"
