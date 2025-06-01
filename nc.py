# -*- coding: utf-8 -*-
import xarray as xr
import pandas as pd
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import numpy as np
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import os

def classify_enso(temp):
    """
    Classify ENSO conditions based on temperature anomaly:
    - Strong El Niño: > 1.5°C
    - Moderate El Niño: 1.0°C to 1.5°C
    - Weak El Niño: 0.5°C to 1.0°C
    - Normal: -0.5°C to 0.5°C
    - Weak La Niña: -1.0°C to -0.5°C
    - Moderate La Niña: -1.5°C to -1.0°C
    - Strong La Niña: < -1.5°C
    """
    if temp > 1.5:
        return "El Niño Kuat"
    elif 1.0 <= temp <= 1.5:
        return "El Niño Sedang"
    elif 0.5 <= temp < 1.0:
        return "El Niño Lemah"
    elif -0.5 <= temp < 0.5:
        return "Normal"
    elif -1.0 <= temp < -0.5:
        return "La Niña Lemah"
    elif -1.5 <= temp < -1.0:
        return "La Niña Sedang"
    else:
        return "La Niña Kuat"

# Paths to both NetCDF files
file_1995_2020 = "data_stream-oper_stepType-instant.nc"
file_2021_2025 = os.path.join(
    "C:/Users/Geral/Downloads/Punya Adin/Surface temp & total precipitation_2021-2025",
    "data_stream-oper_stepType-instant.nc"
)

# Open and combine both datasets
print("Membuka dan menggabungkan data 1995-2020 dan 2021-2025...")
ds1 = xr.open_dataset(file_1995_2020)
ds2 = xr.open_dataset(file_2021_2025)

# Concatenate along the time dimension (valid_time)
ds = xr.concat([ds1, ds2], dim="valid_time")

# Display dataset information
print("\nDataset Information:")
print(ds)

# Display available variables
print("\nAvailable Variables:")
print(ds.variables)

# Display coordinates
print("\nCoordinates:")
print(ds.coords)

# If 'sst' variable exists, show its attributes
if 'sst' in ds:
    print("\nSST Variable Information:")
    print(ds['sst'])

# Get sea surface temperature (sst) data and convert from Kelvin to Celsius
sst = ds['sst'] - 273.15  # Convert from Kelvin to Celsius

# Calculate mean temperature across time dimension
mean_sst = sst.mean(dim='valid_time')

# Create the map
plt.figure(figsize=(15, 10))
ax = plt.axes(projection=ccrs.PlateCarree())

# Set map extent to Pangandaran region with some padding
ax.set_extent([108.2, 108.9, -7.9, -7.5], crs=ccrs.PlateCarree())

# Add detailed features
ax.add_feature(cfeature.LAND, facecolor='lightgray', alpha=0.5)
ax.add_feature(cfeature.OCEAN, facecolor='lightblue', alpha=0.5)
ax.add_feature(cfeature.COASTLINE, linewidth=1)
ax.add_feature(cfeature.RIVERS, linewidth=0.5)
ax.add_feature(cfeature.LAKES, alpha=0.5)

# Add more detailed land features
states_provinces = cfeature.NaturalEarthFeature(
    category='cultural',
    name='admin_1_states_provinces_lines',
    scale='10m',
    facecolor='none')
ax.add_feature(states_provinces, edgecolor='gray', linewidth=0.5)

# Add gridlines with labels
gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                  linewidth=0.5, color='gray', alpha=0.5, linestyle='--')
gl.top_labels = False
gl.right_labels = False
gl.xformatter = LONGITUDE_FORMATTER
gl.yformatter = LATITUDE_FORMATTER
gl.xlabel_style = {'size': 8}
gl.ylabel_style = {'size': 8}

# Plot the temperature data with improved styling
im = mean_sst.plot(ax=ax, transform=ccrs.PlateCarree(),
                  cmap='coolwarm',  # Using a more natural colormap
                  vmin=22, vmax=31,  # Temperature range
                  cbar_kwargs={
                      'label': 'Suhu Permukaan Laut (°C)',
                      'orientation': 'horizontal',
                      'pad': 0.1,
                      'aspect': 40,
                      'shrink': 0.8
                  })

# Add title with better formatting
plt.title('Peta Suhu Permukaan Laut Rata-rata (1995-2025)\nKabupaten Pangandaran, Jawa Barat',
          pad=20, size=14, weight='bold')

# Add a text box with coordinates
plt.text(0.02, 0.02, 'Koordinat: 108.3°E-108.8°E, 7.6°S-7.85°S',
         transform=ax.transAxes, bbox=dict(facecolor='white', alpha=0.7),
         size=8)

# Save the map with high resolution
plt.savefig('peta_suhu_pangandaran_1995_2025.png', dpi=300, bbox_inches='tight')
plt.close()

# Create a DataFrame with the time series data
df = sst.to_dataframe().reset_index()

# Convert valid_time to datetime if it's not already
df['valid_time'] = pd.to_datetime(df['valid_time'])

# Extract year and month
df['year'] = df['valid_time'].dt.year
df['month'] = df['valid_time'].dt.month

# Calculate monthly averages
monthly_avg = df.groupby(['year', 'month'])['sst'].mean().reset_index()

# Calculate the overall mean temperature for the reference period
mean_temp = monthly_avg['sst'].mean()

# Calculate temperature anomalies
monthly_avg['anomaly'] = monthly_avg['sst'] - mean_temp

# Add ENSO classification
monthly_avg['enso'] = monthly_avg['anomaly'].apply(classify_enso)

# Create month names in Indonesian
month_names = {
    1: 'Januari', 2: 'Februari', 3: 'Maret', 4: 'April',
    5: 'Mei', 6: 'Juni', 7: 'Juli', 8: 'Agustus',
    9: 'September', 10: 'Oktober', 11: 'November', 12: 'Desember'
}
monthly_avg['bulan'] = monthly_avg['month'].map(month_names)

# Round temperatures and anomalies to 2 decimal places
monthly_avg['sst'] = monthly_avg['sst'].round(2)
monthly_avg['anomaly'] = monthly_avg['anomaly'].round(2)

# Create Excel writer object
excel_file = 'Suhu_Permukaan_Laut_Bulanan_1995_2025.xlsx'
with pd.ExcelWriter(excel_file, engine='openpyxl') as writer:
    # Create sheets for different views of the data
    
    # Sheet 1: Monthly temperatures
    pivot_temp = monthly_avg.pivot(index='year', columns='bulan', values='sst')
    pivot_temp = pivot_temp[list(month_names.values())]
    pivot_temp.to_excel(writer, sheet_name='Suhu Bulanan', startrow=4)
    
    # Sheet 2: Temperature anomalies
    pivot_anomaly = monthly_avg.pivot(index='year', columns='bulan', values='anomaly')
    pivot_anomaly = pivot_anomaly[list(month_names.values())]
    pivot_anomaly.to_excel(writer, sheet_name='Anomali Suhu', startrow=4)
    
    # Sheet 3: ENSO classification
    pivot_enso = monthly_avg.pivot(index='year', columns='bulan', values='enso')
    pivot_enso = pivot_enso[list(month_names.values())]
    pivot_enso.to_excel(writer, sheet_name='Klasifikasi ENSO', startrow=4)
    
    # Format each sheet
    for sheet_name in ['Suhu Bulanan', 'Anomali Suhu', 'Klasifikasi ENSO']:
        worksheet = writer.sheets[sheet_name]
        
        # Add titles and information
        if sheet_name == 'Suhu Bulanan':
            title = 'Suhu Permukaan Laut Bulanan 1995-2025'
            subtitle = 'Suhu dalam derajat Celsius (°C)'
        elif sheet_name == 'Anomali Suhu':
            title = 'Anomali Suhu Permukaan Laut 1995-2025'
            subtitle = 'Anomali dalam derajat Celsius (°C)'
        else:
            title = 'Klasifikasi ENSO Berdasarkan Anomali Suhu 1995-2025'
            subtitle = 'Klasifikasi: El Niño Kuat/Sedang/Lemah, Normal, La Niña Kuat/Sedang/Lemah'
            
        worksheet.cell(row=1, column=1, value=title)
        worksheet.cell(row=2, column=1, value=subtitle)
        worksheet.cell(row=3, column=1, value='Lokasi: Kabupaten Pangandaran')
        
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

print("\nPeta telah disimpan dalam file: peta_suhu_pangandaran_1995_2025.png")
print("Data telah disimpan dalam file Excel: {}".format(excel_file))
print("File berisi data suhu rata-rata bulanan dan ENSO dari tahun 1995-2025")

# Close the dataset
ds.close()
ds1.close()
ds2.close()
