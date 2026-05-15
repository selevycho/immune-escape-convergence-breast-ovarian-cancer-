import matplotlib.pyplot as plt
import pandas as pd

def main():
    # Final statistics based on our inventory analysis
    # Breast_Cancer: 122
    # Ovarian_Cancer: 96
    # Colon_Cancer: 105
    # Other/Technical: 5
    
    data = {
        'Cancer Type': ['Breast Cancer', 'Ovarian Cancer', 'Colon Cancer', 'Other'],
        'Sample Count': [122, 96, 105, 5]
    }
    
    df = pd.DataFrame(data)
    
    # Define professional colors for the publication
    colors = ['#ff9999', '#66b3ff', '#99ff99', '#ffcc99']
    
    # Create the figure
    plt.figure(figsize=(10, 7))
    
    # Plotting the pie chart
    plt.pie(
        df['Sample Count'], 
        labels=df['Cancer Type'], 
        autopct='%1.1f%%', 
        startangle=140, 
        colors=colors,
        explode=(0.05, 0.05, 0, 0), # Highlight Breast and Ovarian
        shadow=True
    )
    
    plt.title('Distribution of Samples in CPTAC-2 Project (N=328)', fontsize=14, fontweight='bold')
    
    # Add a legend
    plt.legend(df['Cancer Type'], title="Cancer Types", loc="center left", bbox_to_anchor=(1, 0, 0.5, 1))
    
    # Save the plot for the thesis
    plt.tight_layout()
    plt.savefig('results/sample_distribution_pie_chart.png', dpi=300)
    print("Pie chart successfully saved to results/sample_distribution_pie_chart.png")

if __name__ == "__main__":
    main()
