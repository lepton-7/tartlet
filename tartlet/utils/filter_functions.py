from tart.utils.activity_inference import Candidate


def default_check(cand: Candidate):
    issues = []

    if abs(cand.from_switch_end_relative) > 0.2:
        issues.append("Too far from end")
    if cand.rel_cov_delta > -0.2:
        issues.append("Small drop")
    if cand.symks_pval > 0.05:
        issues.append("Symmetry check failed")

    if len(issues) > 0:
        cand.failnote = ". ".join(issues)
        return "fail"

    return "pass"
