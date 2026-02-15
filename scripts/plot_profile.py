import pandas as pd
import matplotlib.pyplot as plt
import sys

def main(csv_file):
    # Load profile data
    df = pd.read_csv(csv_file)
    
    # Sort by Total Peak Memory to visualize the cumulative "weight" of the engine
    #df = df.sort_values(by='TotalPeakMem(KB)', ascending=True)
    
    fig, ax1 = plt.subplots(figsize=(13, 8))

    # 1. Bar chart for Execution Time (Primary Y-axis)
    color_time = 'tab:blue'
    ax1.set_xlabel('Methods (Sorted by Total Memory Footprint)', fontweight='bold')
    ax1.set_ylabel('Time (ms)', color=color_time, fontweight='bold')
    ax1.bar(df['Method'], df['TotalTime(ms)'], color=color_time, alpha=0.5, label='Time (ms)')
    ax1.tick_params(axis='y', labelcolor=color_time)
    
    plt.xticks(rotation=35, ha='right')

    # 2. Line charts for Memory (Secondary Y-axis)
    ax2 = ax1.twinx()
    
    # Differential Memory (The "Spike")
    color_diff = 'tab:red'
    ax2.plot(df['Method'], df['MemDiff(KB)'], color=color_diff, marker='o', 
             linestyle='-', linewidth=2, markersize=8, label='Memory Diff (Delta KB)')
    
    # Absolute Peak Memory (The "Ceiling")
    color_total = 'darkred'
    ax2.plot(df['Method'], df['TotalPeakMem(KB)'], color=color_total, marker='s', 
             linestyle='--', linewidth=2, markersize=6, alpha=0.7, label='Total Peak (RSS KB)')

    ax2.set_ylabel('Memory (KB)', color='tab:red', fontweight='bold')
    ax2.tick_params(axis='y', labelcolor='tab:red')

    # Formatting and Styling
    plt.title('CFEM High-Scale Performance: Time vs. Absolute Memory Scaling', fontsize=15, pad=25)
    
    # Combined Legend
    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(lines1 + lines2, labels1 + labels2, loc='upper left', frameon=True, shadow=True)

    ax1.grid(axis='y', linestyle=':', alpha=0.6)
    
    # Add data labels for the Total Peak to see the GB count clearly
    for i, txt in enumerate(df['TotalPeakMem(KB)']):
        # Only label if it's a significant amount (e.g., > 0.5 GB) to avoid clutter
        if txt > 500000:
            ax2.annotate(f"{txt/1024/1024:.1f} GB", (df['Method'].iloc[i], df['TotalPeakMem(KB)'].iloc[i]),
                         textcoords="offset points", xytext=(0,10), ha='center', fontsize=9, color='darkred')

    fig.tight_layout()
    plt.savefig('performance_report_deep.png', dpi=300)
    print(f"Deep Analysis Report saved as performance_report_deep.png based on {csv_file}")

if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "profile_results.csv")