import pandas as pd
import matplotlib.pyplot as plt
import sys

def main(csv_file):
    # Load profile data
    df = pd.read_csv(csv_file)
    
    # 1. Sort by memory usage (ascending) to visualize scaling trends
    # df = df.sort_values(by='MemDiff(KB)', ascending=True)
    # df = df.sort_values(by='TotalTime(ms)', ascending=True)
    
    fig, ax1 = plt.subplots(figsize=(12, 7))

    # Bar chart for Time (Primary Y-axis)
    color = 'tab:blue'
    ax1.set_xlabel('Methods')
    ax1.set_ylabel('Time (ms)', color=color)
    ax1.bar(df['Method'], df['TotalTime(ms)'], color=color, alpha=0.6, label='Time')
    ax1.tick_params(axis='y', labelcolor=color)
    
    # Rotation adjusted to 30 degrees for longer method names
    plt.xticks(rotation=30, ha='right')

    # Line chart for Memory usage (Secondary Y-axis)
    ax2 = ax1.twinx()
    color = 'tab:red'
    ax2.set_ylabel('Memory Diff (KB)', color=color)
    ax2.plot(df['Method'], df['MemDiff(KB)'], color=color, marker='o', linewidth=2, label='Memory')
    ax2.tick_params(axis='y', labelcolor=color)

    plt.title('FEA Performance Profile: Time vs Memory')
    
    # Add a legend for both axes
    lines, labels = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax2.legend(lines + lines2, labels + labels2, loc='upper left')

    fig.tight_layout()
    plt.savefig('performance_report.png')
    print(f"Report saved as performance_report.png based on {csv_file}")

if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "profile_results.csv")