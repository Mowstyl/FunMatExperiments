import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


exts = ["svg", "png"]
lgnd_loc = (0.03, 0.55)
fig_size = (10.80, 4.93)
labels = ["Deutsch-Jozsa", "_DJ confidence interval", "Bernstein–Vazirani", "_BV confidence interval", "Grover iteration", "_G confidence interval"]

raw = pd.read_csv("raw/dataQJ_1691508388_2_10.csv", sep=';')
raw.rename(columns={"num_qubits": "Number of Qubits",
                    "DJ_times": "Deutsch-Jozsa Execution Time", "DJ_sizes": "Deutsch-Jozsa Memory Usage",
                    "BV_times": "Bernstein–Vazirani Execution Time", "BV_sizes": "Bernstein–Vazirani Memory Usage",
                    "grover_times": "Grover Iteration Execution Time", "grover_sizes": "Grover Iteration Memory Usage"},
           inplace=True)

times = sns.lineplot(data=raw, color='k', legend=True, x="Number of Qubits", y="Deutsch-Jozsa Execution Time")
times = sns.lineplot(data=raw, color='k', legend=True, linestyle="dashed", x="Number of Qubits", y="Bernstein–Vazirani Execution Time", ax=times)
times = sns.lineplot(data=raw, color='k', legend=True, linestyle="dotted", x="Number of Qubits", y="Grover Iteration Execution Time", ax=times)
times.set_ylabel("Time (s)")
times.set_title("Average execution time (full vector)")
times.legend(labels, loc=lgnd_loc)
#h, l = times.get_legend_handles_labels()
#print(h, l)
#times.legend(h[::2], labels[::2], loc=lgnd_loc)

for ext in exts:
    times.get_figure().savefig(f"charts/timesQJ_DJBVG_1-10.{ext}")

plt.clf()

sizes = sns.lineplot(data=raw, color='k', legend=True, x="Number of Qubits", y="Deutsch-Jozsa Memory Usage", ax=None)
sizes = sns.lineplot(data=raw, color='k', legend=True, linestyle="dashed", x="Number of Qubits", y="Bernstein–Vazirani Memory Usage", ax=sizes)
sizes = sns.lineplot(data=raw, color='k', legend=True, linestyle="dotted", x="Number of Qubits", y="Grover Iteration Memory Usage", ax=sizes)
sizes.set_ylabel("Memory (Bytes)")
sizes.set_title("Average max memory usage (full vector)")
sizes.legend(labels, loc=lgnd_loc)
#h, l = sizes.get_legend_handles_labels()
#sizes.legend(h[::2], labels[::2], loc=lgnd_loc)

for ext in exts:
    sizes.get_figure().savefig(f"charts/sizesQJ_DJBVG_1-10.{ext}")


'''
#df = raw.groupby("num_qubits").mean()
#cols = raw.columns[1:4:2]
cols = raw.columns[1:]
for i in range(len(cols)):
    col_name = cols[i]
    ax = sns.lineplot(data=raw, x="Number of Qubits", y=col_name)
    ax.set_ylabel("Time (μs)")
'''
