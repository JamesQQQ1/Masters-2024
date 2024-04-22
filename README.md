# Masters-Research

## Complete Model for Master's Research

### Overview
This repository contains the complete model developed as part of my Master's research. It is designed to run analyses once all data manipulation processes have been completed.

### File Descriptions

- `Prophet.py`: This script utilizes the Prophet forecasting model to generate predictions based on time series data.
- `Raster_Layer.py`: This script handles the conversion of ArcGIS raster files to the NetCDF format.
- `extrapo_population.py`: This script extrapolates population data to estimate population distribution across geographical regions.
- `final_2.6.py`: This script represents one of the final versions of the model, tailored for scenario 2.6.
- `final_4.5.py`: This script represents one of the final versions of the model, tailored for scenario 4.5.
- `final_8.5.py`: This script represents one of the final versions of the model, tailored for scenario 8.5.
- `land_use_change.py`: This script analyzes changes in land use patterns over time.
- `land_use_slice.py`: This script slices and processes land use data to generate inputs for the model.
- `last_year_avg.py`: This script calculates the average values of relevant variables from the last year of data.

### Data Preparation

The model expects input files in NetCDF format. Ensure that all your data files conform to this format or modify the script to accommodate different file formats as needed.

## Model Execution Steps

To run the model for various scenarios, follow these steps in the specified order using the provided scripts:

1. **Prepare Land Use Data:**
   - Run `land_use_slice.py` to slice and process land use data:
     ```
     python land_use_slice.py
     ```
   - Execute `land_use_change.py` to analyze changes in land use patterns over time:
     ```
     python land_use_change.py
     ```
   - Use `last_year_avg.py` to calculate average values of relevant variables from the last year of data:
     ```
     python last_year_avg.py
     ```

2. **Convert Raster Files:**
   - Execute `Raster_Layer.py` to convert ArcGIS raster files to NetCDF format:
     ```
     python Raster_Layer.py
     ```

3. **Forecast Time Series Data:**
   - Run `Prophet.py`, which utilizes the Prophet forecasting model for time series predictions:
     ```
     python Prophet.py
     ```

4. **Extrapolate Population Data:**
   - Use `extrapo_population.py` to estimate population distribution across geographical regions:
     ```
     python extrapo_population.py
     ```

5. **Run Scenario-Specific Models:**
   - Depending on the scenario you are working with, execute the corresponding script:
     - For scenario 2.6:
       ```
       python final_2.6.py
       ```
     - For scenario 4.5:
       ```
       python final_4.5.py
       ```
     - For scenario 8.5:
       ```
       python final_8.5.py
       ```

6. **Validation and Distance Analysis:**
   - Validate the model outputs with `validation.py`:
     ```
     python validation.py
     ```
   - Analyze distance metrics using `distance.py`:
     ```
     python distance.py
     ```

Follow these steps closely to ensure correct preparation and analysis of your data, leading to reliable outputs for each scenario.

### Raw Data Files and Flexibility

#### Included Data Files

This repository includes raw data files that were specifically used in my research. These files serve as examples or starting points for users who wish to understand the data format and structure required by the model.


By adhering to these guidelines, users can effectively leverage the model for a wide array of research questions, making the most of its capabilities to analyze and interpret wind energy potential, land use impacts, or other environmental and geographical phenomena.

For any queries or further assistance, feel free to open an issue in this repository.
