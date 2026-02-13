import pandas as pd
import matplotlib.pyplot as plt
import sys

def main(csv_file):
    # Load profile data
    df = pd.read_csv(csv_file)
    
    # 1. Sort by memory usage (ascending) to visualize scaling trends
    # This highlights the correlation between connectivity storage and time
    df = df.sort_values(by='MemDiff(KB)', ascending=True)
    
    fig, ax1 = plt.subplots(figsize=(12, 7))

    # Bar chart for Time (Primary Y-axis)
    color = 'tab:blue'
    ax1.set_xlabel('Methods (Sorted by Memory Usage)')
    ax1.set_ylabel('Time (ms)', color=color, fontweight='bold')
    bars = ax1.bar(df['Method'], df['TotalTime(ms)'], color=color, alpha=0.6, label='Time (ms)')
    ax1.tick_params(axis='y', labelcolor=color)
    
    # Rotation adjusted for scannability
    plt.xticks(rotation=30, ha='right')

    # Line chart for Memory usage (Secondary Y-axis)
    ax2 = ax1.twinx()
    color = 'tab:red'
    ax2.set_ylabel('Memory Diff (KB)', color=color, fontweight='bold')
    line = ax2.plot(df['Method'], df['MemDiff(KB)'], color=color, marker='o', 
                    linewidth=2, markersize=8, label='Memory (KB)')
    ax2.tick_params(axis='y', labelcolor=color)

    plt.title('CFEM Engine Performance: Time vs Memory Scaling', fontsize=14, pad=20)
    
    # Combined legend for clarity
    lines, labels = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(lines + lines2, labels + labels2, loc='upper left', frameon=True)

    # Grid for background reference
    ax1.grid(axis='y', linestyle=':', alpha=0.5)

    fig.tight_layout()
    plt.savefig('performance_report.png', dpi=300)
    print(f"Report saved as performance_report.png based on {csv_file}")

if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "profile_results.csv")