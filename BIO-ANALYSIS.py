#!/usr/bin/env python
# coding: utf-8

# Ανάλυση Δυναμικής MRI για Ταξινόμηση Συμπεριφοράς Εικονοστοιχείων και Οπτικοποίηση (NIfTI & GIF)

# In[ ]:


import SimpleITK as sitk
import numpy as np
import matplotlib.pyplot as plt
import os
import imageio 
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch
import traceback 

MRI_DIR = r"E:\Biosignals\mris"
MASK_DIR = r"E:\Biosignals\annotations" 
OUTPUT_DIR = r"E:\Biosignals\Nifti" 

label_dict = {
    1: "Uptake",    
    2: "Plateau",   
    3: "Washout"    
}

color_dict = {
    1: (1, 0, 0),    
    2: (0, 1, 0),    
    3: (0, 0, 1)     
}

colors_for_nifti_map = [
    (0, 0, 0, 0),
    (1, 0, 0, 1),
    (0, 1, 0, 1),
    (0, 0, 1, 1)
]
nifti_cmap = ListedColormap(colors_for_nifti_map)


def classify_pixel_behavior(image0000_np, image0001_np, roi_mask_np, threshold_percent=10):
    behavior_map_np = np.zeros_like(image0000_np, dtype=np.uint8)
    roi_indices = np.where(roi_mask_np > 0)

    for z, y, x in zip(*roi_indices):
        intensity_t0 = float(image0000_np[z, y, x])
        intensity_t1 = float(image0001_np[z, y, x])

        if intensity_t0 < 1e-6: 
            if intensity_t1 > 1e-6:
                behavior_map_np[z, y, x] = 1 
            else:
                behavior_map_np[z, y, x] = 2 
        else:
            percentage_change = ((intensity_t1 - intensity_t0) / intensity_t0) * 100

            if percentage_change > threshold_percent:
                behavior_map_np[z, y, x] = 1  
            elif -threshold_percent <= percentage_change <= threshold_percent:
                behavior_map_np[z, y, x] = 2  
            elif percentage_change < -threshold_percent:
                behavior_map_np[z, y, x] = 3  

    return behavior_map_np


def create_overlay_gif(image1, classified_behavior_map, output_path, current_label_value, label_name, specific_rgb_color):
    frames = []
    mask_for_gif = (classified_behavior_map == current_label_value).astype(np.uint8)
    overlay_alpha = 0.6

    custom_cmap_colors = [
        (0, 0, 0, 0),  
        (*specific_rgb_color, overlay_alpha)
    ]
    custom_cmap = ListedColormap(custom_cmap_colors)

    for z in range(image1.shape[0]):
        if np.sum(mask_for_gif[z]) == 0:
            continue

        fig, ax = plt.subplots(figsize=(6, 6))
        ax.imshow(image1[z], cmap='gray')
        ax.imshow(mask_for_gif[z], cmap=custom_cmap, vmin=0, vmax=1)
        ax.axis('off') 
        ax.set_title(f"{label_name} - Slice {z}", fontsize=12)

        legend_color_rgba = (*specific_rgb_color, 0.7) 
        legend_elements = [Patch(facecolor=legend_color_rgba, edgecolor='none', label=label_name)]
        ax.legend(handles=legend_elements, loc='lower right', fontsize=10, framealpha=0.7)

        fig.canvas.draw()
        frame = np.array(fig.canvas.buffer_rgba())
        frame = frame.reshape(fig.canvas.get_width_height()[::-1] + (4,)) 
        frame_rgb = frame[:, :, :3] 
        frames.append(frame_rgb)
        plt.close(fig) 

    if frames: 
        imageio.mimsave(output_path, frames, duration=0.2)
        print(f"✅ Saved GIF to: {output_path}")
    else:
        print(f"Skipped GIF creation for {label_name}: No relevant regions found.")


def process_patient(mri_folder, annot_folder, output_folder, label_map, color_map_for_gif):
    os.makedirs(output_folder, exist_ok=True)
    output_nifti_folder = os.path.join(output_folder, "nifti_behavior_maps")
    os.makedirs(output_nifti_folder, exist_ok=True)
    output_gifs_folder = os.path.join(output_folder, "gifs")
    os.makedirs(output_gifs_folder, exist_ok=True)

    for file in os.listdir(mri_folder):
        if file.endswith("_0000.nii.gz"): 
            patient_id = file.replace("_0000.nii.gz", "")
            print(f"\n🧪 Processing patient: {patient_id}")
            try:
                image0000_path = os.path.join(mri_folder, f"{patient_id}_0000.nii.gz")
                image0001_path = os.path.join(mri_folder, f"{patient_id}_0001.nii.gz")
                roi_mask_path = os.path.join(annot_folder, f"{patient_id}.nii.gz") 

                image0000_sitk = sitk.ReadImage(image0000_path)
                image0001_sitk = sitk.ReadImage(image0001_path)
                roi_mask_sitk = sitk.ReadImage(roi_mask_path)

                image0000_np = sitk.GetArrayFromImage(image0000_sitk)
                image0001_np = sitk.GetArrayFromImage(image0001_sitk)
                roi_mask_np = sitk.GetArrayFromImage(roi_mask_sitk)

                if not (image0000_np.shape == image0001_np.shape == roi_mask_np.shape):
                    print(f"  ❌ Dimension mismatch for {patient_id}. Skipping.")
                    continue

                print("  Classifying pixel behavior...")
                classified_behavior_map_np = classify_pixel_behavior(
                    image0000_np, image0001_np, roi_mask_np, threshold_percent=10
                )
                print("  Classification complete.")

                               
                patient_nifti_folder = os.path.join(output_nifti_folder, patient_id)
                os.makedirs(patient_nifti_folder, exist_ok=True)
                
                nifti_filename = "behavior_map.nii.gz"
                nifti_full_path = os.path.join(patient_nifti_folder, nifti_filename)

                
                behavior_map_sitk = sitk.GetImageFromArray(classified_behavior_map_np)
                behavior_map_sitk.CopyInformation(image0000_sitk) 
                sitk.WriteImage(behavior_map_sitk, nifti_full_path)
                print(f"✅ Saved NIfTI to: {nifti_full_path}")

                print("  Generating GIFs...")
                for label_value, label_name in label_map.items():
                    gif_filename = f"{patient_id}_{label_name}.gif"
                    gif_full_path = os.path.join(output_gifs_folder, gif_filename)
                    specific_rgb_color = color_map_for_gif.get(label_value, (0.5, 0.5, 0.5))

                    create_overlay_gif(
                        image0000_np, 
                        classified_behavior_map_np, 
                        gif_full_path, 
                        label_value, 
                        label_name, 
                        specific_rgb_color
                    )
                print("  GIFs complete.")

            except Exception as e:
                print(f"  ❌ Error for {patient_id}: {e}")
                traceback.print_exc() 


if __name__ == "__main__":
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    os.makedirs(os.path.join(OUTPUT_DIR, "nifti_behavior_maps"), exist_ok=True)
    os.makedirs(os.path.join(OUTPUT_DIR, "gifs"), exist_ok=True)

    print("Starting analysis...")
    process_patient(MRI_DIR, MASK_DIR, OUTPUT_DIR, label_dict, color_dict)
    print("\n--- Analysis finished ---")
    print(f"Results saved in: {OUTPUT_DIR}")


# Εξαγωγή και Οπτικοποίηση Χαρακτηριστικών Ακτινομικής (Radiomics) για Ανάλυση Ιατρικών Εικόνων

# In[44]:


import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import SimpleITK as sitk
from radiomics import featureextractor
from sklearn import manifold
from sklearn.metrics import euclidean_distances
from tqdm import tqdm

# Set your paths here
MRI_DIR = r"E:\Biosignals\mris"
MASK_DIR = r"E:\Biosignals\annotations"
OUTPUT_DIR = r"E:\Biosignals\output"
os.makedirs(OUTPUT_DIR, exist_ok=True)

# 1. File Handling and Feature Extraction
def extract_features(mri_dir, mask_dir):
    # Initialize feature extractor with custom settings
    settings = {
        'binWidth': 25,
        'interpolator': 'sitkBSpline',
        'resampledPixelSpacing': [1, 1, 1],
        'geometryTolerance': 1e-6,
        'normalize': False
    }
    extractor = featureextractor.RadiomicsFeatureExtractor(**settings)
    
    # Get all case files
    mri_files = sorted(glob.glob(os.path.join(mri_dir, '*.nii.gz')))
    mask_files = sorted(glob.glob(os.path.join(mask_dir, '*.nii.gz')))
    
    features = {}
    for mri_path, mask_path in tqdm(zip(mri_files, mask_files), desc="Extracting features"):
        try:
            case_id = os.path.basename(mri_path).split('.')[0]
            image = sitk.ReadImage(mri_path)
            mask = sitk.ReadImage(mask_path)
            
            # Ensure mask is binary (0 and 1)
            mask = sitk.BinaryThreshold(mask, lowerThreshold=1, upperThreshold=255)
            
            features[case_id] = extractor.execute(image, mask)
        except Exception as e:
            print(f"Error processing {case_id}: {str(e)}")
    
    return features

# 2. Feature Processing
def prepare_feature_matrix(features):
    # Get valid feature names (excluding diagnostics)
    feature_names = sorted([k for k in features[next(iter(features))].keys() 
                          if k.startswith("original_")])
    
    # Create numpy array of all values
    samples = np.zeros((len(features), len(feature_names)))
    for i, case_id in enumerate(features.keys()):
        samples[i,:] = [features[case_id][fn] for fn in feature_names]
    
    # Handle NaNs if any
    samples = np.nan_to_num(samples)
    
    return samples, feature_names

# 3. Visualization Functions
def plot_mds(samples, feature_names):
    """Enhanced MDS plot matching the first code's style and groups"""
    # Calculate similarity matrix
    similarities = euclidean_distances(samples)
    
    # Apply MDS with the same parameters
    mds = manifold.MDS(n_components=2, max_iter=5000, eps=1e-12,
                      dissimilarity="precomputed", n_jobs=1, metric=False,
                      random_state=3)  # Added random_state for reproducibility
    pos = mds.fit_transform(similarities)
    
    # Set up the plot with the same style
    fig = plt.figure(1, figsize=(10, 8))
    ax = plt.axes([0., 0., 1., 1.])
    s = 100  # Marker size
    
    # Define the same groups as in the first code
    meningioma = [0, 1, 2]
    glioma = [3, 5, 9]
    astrocytoma = [4, 6, 7, 8]
    
    # Plot with the same colors and styling
    plt.scatter(pos[meningioma, 0], pos[meningioma, 1], 
               color='navy', alpha=1.0, s=s, lw=1, label='meningioma')
    plt.scatter(pos[glioma, 0], pos[glioma, 1], 
               color='turquoise', alpha=1.0, s=s, lw=1, label='glioma')
    plt.scatter(pos[astrocytoma, 0], pos[astrocytoma, 1], 
               color='darkorange', alpha=0.5, s=s, lw=1, label='astrocytoma')
    
    # Add legend in the same position
    plt.legend(scatterpoints=1, loc=5, shadow=False)
    
    # Calculate and normalize similarities like in the first code
    with np.errstate(divide='ignore'):  # Ignore division by zero warnings
        similarities_normalized = similarities.max() / similarities * 100
        similarities_normalized[np.isinf(similarities_normalized)] = 0
    
    plt.title("MDS Visualization of Radiomics Features")
    plt.xlabel("Component 1")
    plt.ylabel("Component 2")
    
    # Save and show
    plt.savefig(os.path.join(OUTPUT_DIR, 'mds_plot.png'))
    plt.show()
    
    return similarities_normalized  # Optionally return the similarities matrix

def plot_feature_heatmap(df):
    corr = df.corr()
    plt.figure(figsize=(15, 10))
    sns.heatmap(corr, vmax=.8, square=True)
    plt.title("Feature Correlation Heatmap")
    plt.savefig(os.path.join(OUTPUT_DIR, 'feature_heatmap.png'))
    plt.show()

def plot_clustermap(df):
    plt.figure(figsize=(13, 13))
    pp = sns.clustermap(df.corr(), linewidths=.5, figsize=(13,13))
    plt.setp(pp.ax_heatmap.get_yticklabels(), rotation=0)
    plt.savefig(os.path.join(OUTPUT_DIR, 'feature_clustermap.png'))
    plt.show()

# Main Analysis Pipeline
def run_analysis():
    print("Starting radiomics analysis...")
    
    # 1. Extract features
    features = extract_features(MRI_DIR, MASK_DIR)
    
    # 2. Prepare feature matrix
    samples, feature_names = prepare_feature_matrix(features)
    df = pd.DataFrame(samples, columns=feature_names)
    
    # 3. Save raw features
    df.to_csv(os.path.join(OUTPUT_DIR, 'radiomics_features.csv'))
    
    # 4. Visualizations
    print("Creating visualizations...")
    plot_mds(samples, feature_names)
    plot_feature_heatmap(df)
    
    # For clustermap, use a subset of features if you have many
    if len(feature_names) > 50:
        plot_clustermap(df.iloc[:, :50])  # First 50 features
    else:
        plot_clustermap(df)
    
    print(f"Analysis complete! Results saved in {OUTPUT_DIR}")

if __name__ == "__main__":
    run_analysis()


# Ομογενοποίηση σημάτων από πολυκεντρικά δεδομένα χρησιμοποιώντας την τεχνική combat

# In[ ]:


import os
import numpy as np
import pandas as pd
import SimpleITK as sitk
from radiomics import featureextractor
from neuroCombat import neuroCombat
from tqdm import tqdm
import glob

MRI_DIR = r"E:\Biosignals\mris"
MASK_DIR = r"E:\Biosignals\annotations"
OUTPUT_DIR = r"E:\Biosignals\output"
os.makedirs(OUTPUT_DIR, exist_ok=True)

def extract_features(mri_dir, mask_dir):
    settings = {
        'binWidth': 25,
        'interpolator': 'sitkBSpline',
        'resampledPixelSpacing': [1, 1, 1],
        'geometryTolerance': 1e-6,
        'normalize': False
    }
    extractor = featureextractor.RadiomicsFeatureExtractor(**settings)
    mri_files = sorted(glob.glob(os.path.join(mri_dir, '*.nii.gz')))
    mask_files = sorted(glob.glob(os.path.join(mask_dir, '*.nii.gz')))
    mri_basenames = {os.path.basename(f).replace('.nii.gz', ''): f for f in mri_files}
    mask_basenames = {os.path.basename(f).replace('.nii.gz', ''): f for f in mask_files}
    common_ids = sorted(list(mri_basenames.keys() & mask_basenames.keys()))
    mri_files_matched = [mri_basenames[uid] for uid in common_ids]
    mask_files_matched = [mask_basenames[uid] for uid in common_ids]
    features = {}
    for mri_path, mask_path in tqdm(zip(mri_files_matched, mask_files_matched), total=len(mri_files_matched)):
        try:
            case_id = os.path.basename(mri_path).replace('.nii.gz', '')
            feature_vector = extractor.execute(mri_path, mask_path)
            filtered_features = {k: v for k, v in feature_vector.items() if k.startswith('original_')}
            features[case_id] = filtered_features
        except:
            continue
    return features

def prepare_feature_matrix(features_dict):
    df_features = pd.DataFrame.from_dict(features_dict, orient='index')
    if df_features.isnull().sum().sum() > 0:
        df_features = df_features.fillna(df_features.mean())
    feature_names = df_features.columns.tolist()
    data = df_features.values.T
    return data, feature_names

extracted_features_dict = extract_features(MRI_DIR, MASK_DIR)
data, feature_names = prepare_feature_matrix(extracted_features_dict)

feature_variances = np.var(data, axis=1)
constant_features_indices = np.where(feature_variances < np.finfo(float).eps)[0]
if len(constant_features_indices) > 0:
    data = np.delete(data, constant_features_indices, axis=0)
    feature_names = [feature_names[i] for i in range(len(feature_names)) if i not in constant_features_indices]

num_samples = data.shape[1]
if num_samples >= 40 and num_samples % 4 == 0:
    samples_per_batch = num_samples // 4
    batch_values = []
    for i in range(1, 5):
        batch_values.extend([i] * samples_per_batch)
elif num_samples >= 2:
    batch_values = [(i % 2) + 1 for i in range(num_samples)]
else:
    raise ValueError("Not enough samples.")

gender_values = [1] * num_samples
covars = pd.DataFrame({'batch': batch_values, 'gender': gender_values})
batch_col = 'batch'
categorical_cols = []

if len(covars[batch_col].unique()) < 2:
    exit()

data_combat = neuroCombat(dat=data, covars=covars, batch_col=batch_col, categorical_cols=categorical_cols)["data"]
output_harmonized_path = os.path.join(OUTPUT_DIR, "harmonized_biosignals_features.csv")
df_harmonized = pd.DataFrame(data_combat.T, columns=feature_names)
df_harmonized.to_csv(output_harmonized_path, index=False)

