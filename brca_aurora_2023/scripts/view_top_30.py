import pandas as pd
import os

# Define paths to your newly created tables
RESULTS_DIR = os.path.expanduser("~/immune_escape_project/brca_aurora_2023/results")
high_file = os.path.join(RESULTS_DIR, "high_silencing_aurora_all.csv")
low_file = os.path.join(RESULTS_DIR, "low_silencing_aurora_all.csv")

def view_data():
    # Check if files exist
    if not os.path.exists(high_file) or not os.path.exists(low_file):
        print("❌ Error: CSV files not found. Make sure you ran the split tables script!")
        return

    # Load data
    high_df = pd.read_csv(high_file)
    low_df = pd.read_csv(low_file)

    # Function to print tables nicely
    def print_section(title, df):
        print("\n" + "="*70)
        print(f"🚀 {title} (Showing TOP 30 by Mutation Count)")
        print("="*70)
        # Use to_string to ensure all columns are visible in the terminal
        print(df.head(30).to_string(index=False))
        print("-" * 70)
        print(f"Total genes in this category: {len(df)}")

    # Display High Silencing
    print_section("HIGH SILENCING (STEALTH GENES)", high_df)
    
    # Display Low Silencing
    print_section("LOW SILENCING (EXPOSED GENES)", low_df)

if __name__ == "__main__":
    view_data()
