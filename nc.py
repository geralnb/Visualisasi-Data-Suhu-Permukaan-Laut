# -*- coding: utf-8 -*-
import xarray as xr
import pandas as pd
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import numpy as np
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import os
from datetime import datetime

# Create output directory for maps
os.makedirs('peta_suhu_jawa', exist_ok=True)

# Open NetCDF file for Java
ds_jawa = xr.open_dataset("pulau-jawa.nc")

# Get sea surface temperature (sst) data and convert from Kelvin to Celsius
sst = ds_jawa['sst'] - 273.15  # Convert from Kelvin to Celsius

# Create month names in Indonesian
month_names = {
    1: 'Januari', 2: 'Februari', 3: 'Maret', 4: 'April',
    5: 'Mei', 6: 'Juni', 7: 'Juli', 8: 'Agustus',
    9: 'September', 10: 'Oktober', 11: 'November', 12: 'Desember'
}

def create_monthly_map(ax, data, month_name):
    """Create a map for the given data on the specified axis"""
    ax.set_extent([105, 115, -8.5, -5.5], crs=ccrs.PlateCarree())
    
    # Add detailed features
    ax.add_feature(cfeature.LAND, facecolor='lightgray', alpha=0.5)
    ax.add_feature(cfeature.OCEAN, facecolor='lightblue', alpha=0.5)
    ax.add_feature(cfeature.COASTLINE, linewidth=0.5)
    ax.add_feature(cfeature.RIVERS, linewidth=0.3)
    ax.add_feature(cfeature.LAKES, alpha=0.5)
    
    # Add more detailed land features
    states_provinces = cfeature.NaturalEarthFeature(
        category='cultural',
        name='admin_1_states_provinces_lines',
        scale='10m',
        facecolor='none')
    ax.add_feature(states_provinces, edgecolor='gray', linewidth=0.3)
    
    # Add gridlines with labels
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                     linewidth=0.3, color='gray', alpha=0.5, linestyle='--')
    gl.top_labels = False
    gl.right_labels = False
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    gl.xlabel_style = {'size': 6}
    gl.ylabel_style = {'size': 6}
    
    # Plot the temperature data
    im = data.plot(ax=ax, transform=ccrs.PlateCarree(),
                  cmap='coolwarm',
                  vmin=22, vmax=31,
                  add_colorbar=False)
    
    # Add month name as title
    ax.set_title(month_name, size=8, pad=2)
    
    return im

# Process data year by year
years = pd.DatetimeIndex(sst.valid_time.values).year.unique()

for year in years:
    print(f"Processing year {year}...")
    
    # Create figure with 12 subplots (4x3 grid)
    fig = plt.figure(figsize=(20, 25))
    
    # Add main title
    fig.suptitle(f'Peta Suhu Permukaan Laut Pulau Jawa Tahun {year}', 
                 fontsize=16, weight='bold', y=0.95)
    
    # Create subplots for each month
    last_im = None  # Store the last image for colorbar
    
    for month in range(1, 13):
        # Get data for this month
        month_data = sst.sel(valid_time=f"{year}-{month:02d}")
        
        if len(month_data.valid_time) > 0:  # Check if we have data for this month
            # Create subplot
            ax = plt.subplot(4, 3, month, projection=ccrs.PlateCarree())
            
            # Create map
            im = create_monthly_map(ax, month_data.mean(dim='valid_time'), month_names[month])
            last_im = im
    
    # Add colorbar at the bottom
    if last_im is not None:
        cbar_ax = fig.add_axes([0.2, 0.05, 0.6, 0.02])
        cbar = plt.colorbar(last_im, cax=cbar_ax, orientation='horizontal',
                          label='Suhu Permukaan Laut (°C)')
    
    # Add text box with coordinates
    fig.text(0.02, 0.02, 'Pulau Jawa: 105°E-115°E, 5.5°S-8.5°S',
             bbox=dict(facecolor='white', alpha=0.7), size=8)
    
    # Adjust layout
    plt.tight_layout(rect=[0, 0.08, 1, 0.95])
    
    # Save the figure
    output_file = os.path.join('peta_suhu_jawa', f'peta_suhu_jawa_{year}.png')
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()

# Create a summary Excel file with monthly averages
print("\nCreating Excel summary...")

# Create a DataFrame with the time series data
df = sst.to_dataframe().reset_index()

# Convert valid_time to datetime if it's not already
df['valid_time'] = pd.to_datetime(df['valid_time'])

# Extract year and month
df['year'] = df['valid_time'].dt.year
df['month'] = df['valid_time'].dt.month

# Calculate monthly averages
monthly_avg = df.groupby(['year', 'month'])['sst'].mean().reset_index()
monthly_avg['bulan'] = monthly_avg['month'].map(month_names)

# Round temperatures to 2 decimal places
monthly_avg['sst'] = monthly_avg['sst'].round(2)

# Create Excel writer object
excel_file = 'Suhu_Permukaan_Laut_Jawa_Bulanan.xlsx'
with pd.ExcelWriter(excel_file, engine='openpyxl') as writer:
    # Create pivot table
    pivot_table = monthly_avg.pivot(index='year', columns='bulan', values='sst')
    pivot_table = pivot_table[list(month_names.values())]
    
    # Write to Excel
    pivot_table.to_excel(writer, sheet_name='Suhu Bulanan', startrow=4)
    
    # Get the worksheet
    worksheet = writer.sheets['Suhu Bulanan']
    
    # Add title and information
    worksheet.cell(row=1, column=1, value='Suhu Permukaan Laut Bulanan Pulau Jawa')
    worksheet.cell(row=2, column=1, value='Suhu dalam derajat Celsius (°C)')
    worksheet.cell(row=3, column=1, value='Lokasi: Pulau Jawa (105°E-115°E, 5.5°S-8.5°S)')
    
    # Adjust column widths
    for column in worksheet.columns:
        max_length = 0
        column = [cell for cell in column]
        for cell in column:
            try:
                if len(str(cell.value)) > max_length:
                    max_length = len(str(cell.value))
            except:
                pass
        adjusted_width = (max_length + 2)
        worksheet.column_dimensions[column[0].column_letter].width = adjusted_width

print(f"\nPeta tahunan telah disimpan dalam folder: peta_suhu_jawa")
print(f"Data rata-rata bulanan telah disimpan dalam file: {excel_file}")

# Close the dataset
ds_jawa.close()
