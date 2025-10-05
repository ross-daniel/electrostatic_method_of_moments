from collections import defaultdict

class ErrorTableCollector:
    def __init__(self):
        self.data = defaultdict(list)  # {"PM": [(N, w_over_h, err), ...], "Galerkin": [...]}

    def add(self, method: str, w_over_h, N: int, error):
        method = method.strip()
        assert method in ("PM", "Galerkin"), f"Unknown method {method}"
        self.data[method].append((int(N), float(w_over_h), float(error)))

    def _sorted(self, method):
        return sorted(self.data.get(method, []), key=lambda t: (t[0], t[1]))

    def to_latex(self, method: str, caption_suffix: str = "errors"):
        rows = self._sorted(method)
        lines = []
        lines.append(r"\begin{table}[h]")
        lines.append(r"\centering")
        lines.append(r"\caption{%s %s}" % (method, caption_suffix))
        lines.append(r"\begin{tabular}{r r r}")
        lines.append(r"\hline")
        lines.append(r"$N$ & $w/h$ & error \\")
        lines.append(r"\hline")
        for N, wh, e in rows:
            wh_str = f"{wh:.6g}"
            e_str  = f"{e:.3e}" if (abs(e) > 1e3 or (e != 0 and abs(e) < 1e-3)) else f"{e:.6g}"
            lines.append(f"{N} & {wh_str} & {e_str} \\\\")
        lines.append(r"\hline")
        lines.append(r"\end{tabular}")
        lines.append(r"\label{tab:%s-errors}" % method.lower())
        lines.append(r"\end{table}")
        return "\n".join(lines)

    def save_csvs(self, base: str = "."):
        for method in ("PM", "Galerkin"):
            rows = self._sorted(method)
            if rows:
                with open(f"{base}/{method.lower()}_errors.csv", "w") as f:
                    f.write("N,w_over_h,error\n")
                    for N, wh, e in rows:
                        f.write(f"{N},{wh},{e}\n")

# ---- Example integration ----
#collector = ErrorTableCollector()

# Wherever you currently print, do this instead:
# print(f"{method} Error for w/h = {wh}, N = {N}: [{err}]")
#def record_error(method, wh, N, err):
    #collector.add(method, wh, N, err)
    # you can keep the print if you still want console output
    # print(f"{method} Error for w/h = {wh}, N = {N}: [{err}]")

# Your nested loops:
# for N in Ns:
#     for wh in w_over_h_values:
#         err_pm = ...
#         record_error("PM", wh, N, err_pm)
#         err_gk = ...
#         record_error("Galerkin", wh, N, err_gk)

# After the loops:
#print("% PM table:\n")
#print(collector.to_latex("PM"))
#print("\n\n% Galerkin table:\n")
#print(collector.to_latex("Galerkin"))

# Optional CSV export:
#collector.save_csvs()