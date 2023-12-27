# DNA Clustering Application
This is a DNA clustering application built using Python and the Tkinter library for the graphical user interface. The application allows users to perform clustering on a set of DNA sequences using various clustering algorithms, visualize the results, and analyze the clusters. Furthermore, the codebase also includes datasets sourced from Kaggle to enhance its practical applications. You can find the Kaggle [link1](https://www.kaggle.com/datasets/neelvasani/humandnadata), [link2](https://www.kaggle.com/datasets/neelvasani/chimpanzee-and-dog-dna) as provided here.

## Table of Contents

- [Installation](#installation)
- [Usage](#usage)
- [Features](#features)
- [Dependencies](#dependencies)

## Installation

To use the DNA Clustering Application, follow these steps:

1. **Clone the Repository:**
   ```bash
   git clone https://github.com/annupriy/DNA-Sequences-Clustering.git

2. **Run the dna_version4.py file.**

## Usage

Follow these steps to use the DNA Clustering Application:

1. **Launch the Application:**
   - Execute the application using the provided script.
   
2. **Load DNA Sequences:**
   - Click on the "Select DNA File" button to load a file containing DNA sequences.
   
3. **Choose Clustering Algorithm:**
   - Select a clustering algorithm from the available options.
   - Set the relevant parameters for the chosen algorithm.

4. **Run Clustering:**
   - Click the corresponding "Run Clustering" button to start the clustering process.

5. **View Results:**
   - Explore the clustering results, silhouette score, and other information presented in the application window.

6. **Explore Additional Functionalities:**
   - Utilize additional functionalities such as:
     - Plotting clusters to visualize the distribution of DNA sequences.
     - Displaying dendrograms for hierarchical clustering results.
     - and more.

Feel free to experiment with different algorithms, adjust parameters, and explore the various features offered by the application.

For detailed instructions on each feature, refer to the respective sections in the application or consult the [Features](#features) section in this README.

## Features

- **KMeans Clustering:** Perform clustering using the KMeans algorithm with a specified number of clusters.
- **KMedoids Clustering:** Utilize the KMedoids algorithm for clustering with a user-defined cluster count.
- **DBSCAN Clustering:** Apply DBSCAN clustering with adjustable parameters for minimum points and neighborhood radius.
- **Agglomerative Clustering:** Perform hierarchical agglomerative clustering with a specified number of clusters.
- **Dendrogram Display:** Visualize hierarchical clustering results with a dendrogram.
- **Cluster Plotting:** Generate plots to visualize the distribution of DNA sequences across clusters.
- **Silhouette Score Calculation:** Evaluate the quality of clustering results using silhouette scores.

## Dependencies

- `tkinter`: Python's standard GUI library.
- `PIL`: Python Imaging Library, used for image-related operations.
- `numpy`: Provides support for working with arrays and matrices.
- `scikit-learn`: A machine learning library for various clustering algorithms.
- `Levenshtein`: A Python extension for computing Levenshtein distances.
- `matplotlib`: A plotting library for creating visualizations.
- `inflect`: A library for converting numbers into words.
