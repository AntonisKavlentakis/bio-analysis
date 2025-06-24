MRI Analysis Pipeline – Dynamic Behavior Classification, Radiomics Feature Extraction & Combat Harmonization
Περιγραφή
Το παρόν pipeline εκτελεί μια ολοκληρωμένη ανάλυση δυναμικών MRI (μορφή .nii.gz) με στόχους:

Ταξινόμηση της συμπεριφοράς voxel (Uptake, Plateau, Washout) σε δύο χρονικές στιγμές (t0, t1) και αποθήκευση σε μορφή NIfTI.

Οπτικοποίηση ταξινομημένων περιοχών με κινούμενες εικόνες τύπου GIF (ανά κατηγορία συμπεριφοράς).

Εξαγωγή χαρακτηριστικών Radiomics από τις περιοχές ενδιαφέροντος (ROIs) για κάθε εικόνα.

Οπτικοποίηση Radiomics δεδομένων με MDS, heatmaps και clustermaps.

Ομογενοποίηση Radiomics χαρακτηριστικών με τη μέθοδο Combat για προσαρμογή batch effects μεταξύ διαφορετικών κέντρων.

E:\Biosignals\
├── mris\              # Περιέχει εικόνες MRI (π.χ. patient_0000.nii.gz, patient_0001.nii.gz)
├── annotations\       # Περιέχει μάσκες ROI (π.χ. patient.nii.gz)
├── Nifti\             # Έξοδοι από ταξινόμηση και GIFs
├── output\            # Έξοδοι από radiomics και combat analysis



Μέρος 1: Ταξινόμηση Συμπεριφοράς Voxels & Οπτικοποίηση
Περιγραφή
Για κάθε voxel εντός ROI, η ταξινόμηση γίνεται ως εξής:

Uptake: σημαντική αύξηση έντασης > 10%

Plateau: μεταβολή μέσα στο ±10%

Washout: μείωση > 10%

Αποθηκεύονται:

NIfTI χάρτης συμπεριφοράς (behavior_map.nii.gz)

GIF εικόνες ανά κατηγορία (Uptake, Plateau, Washout)

Εκτέλεση
Το script process_patient() επεξεργάζεται όλους τους φακέλους:
if __name__ == "__main__":
    process_patient(MRI_DIR, MASK_DIR, OUTPUT_DIR, label_dict, color_dict)
Μέρος 2: Radiomics Feature Extraction & Visualization
Περιγραφή
Εξαγωγή χαρακτηριστικών με Pyradiomics από κάθε MRI + μάσκα.

Αποθήκευση των χαρακτηριστικών σε .csv.

Οπτικοποίηση με:

MDS (Multidimensional Scaling) για προβολή σε 2D χώρο

Heatmap συσχετίσεων χαρακτηριστικών

Clustermap για εύρεση δομών στις συσχετίσεις

Εκτέλεση
python
Αντιγραφή
Επεξεργασία
if __name__ == "__main__":
    run_analysis()
Μέρος 3: Combat Harmonization
Περιγραφή
Η τεχνική Combat χρησιμοποιείται για να εξουδετερώσει τις διαφορές μεταξύ δεδομένων από διαφορετικά "batch" (π.χ. νοσοκομεία).

Βήματα:
Εξαγωγή Radiomics χαρακτηριστικών όπως παραπάνω

Συλλογή metadata (batch labels, π.χ. νοσοκομείο)

Εφαρμογή neuroCombat

python
Αντιγραφή
Επεξεργασία
from neuroCombat import neuroCombat
Η ομογενοποιημένη έξοδος αποθηκεύεται για περαιτέρω ανάλυση ή ταξινόμηση.

Εξαρτήσεις
Πριν την εκτέλεση, εγκατάστησε τις ακόλουθες βιβλιοθήκες:

bash
Αντιγραφή
Επεξεργασία
pip install SimpleITK matplotlib imageio seaborn scikit-learn pandas tqdm pyradiomics neuroCombat
Ρυθμίσεις Pyradiomics
python
Αντιγραφή
Επεξεργασία
settings = {
    'binWidth': 25,
    'interpolator': 'sitkBSpline',
    'resampledPixelSpacing': [1, 1, 1],
    'geometryTolerance': 1e-6,
    'normalize': False
}
Παραγόμενα Αρχεία
Nifti/

nifti_behavior_maps/patient_x/behavior_map.nii.gz

gifs/patient_x_Uptake.gif, ...Plateau.gif, ...Washout.gif

output/

radiomics_features.csv

mds_plot.png

feature_heatmap.png

feature_clustermap.png

Ομάδες Όγκων για Visualization
Χρησιμοποιούνται ενδεικτικές ομάδες για το MDS plot:

meningioma = [0, 1, 2]

glioma = [3, 5, 9]

astrocytoma = [4, 6, 7, 8]

Μπορείς να τις αλλάξεις ανάλογα με τα δεδομένα σου.
