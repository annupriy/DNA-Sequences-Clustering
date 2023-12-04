import tkinter as tk
from tkinter import filedialog
from PIL import Image, ImageTk
import numpy as np
from sklearn.cluster import KMeans
from sklearn_extra.cluster import KMedoids
from sklearn.cluster import DBSCAN
from sklearn.model_selection import train_test_split
import Levenshtein as lev
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from sklearn.metrics import silhouette_score
import inflect
from sklearn.cluster import AgglomerativeClustering
from scipy.cluster import hierarchy
from threading import Thread

class DNAClusteringApp:
    def __init__(self, root):
        self.root = root
        self.root.title("DNA Clustering")

        # Set the background image
        bg_image_path = r"DNA-sequences-clustering-main\DNA-sequences-clustering-main\bg.png"  
        self.bg_image = ImageTk.PhotoImage(Image.open(bg_image_path))
        self.background_label = tk.Label(self.root, image=self.bg_image)
        self.background_label.place(relwidth=1, relheight=1)

        self.dna_sequences = []
        self.true_labels = []
        self.cluster_labels = []

        # Create GUI components
        self.create_widgets()

    def create_widgets(self):
        # File selection button for DNA sequences
        self.select_file_button = tk.Button(self.root, text="Select DNA File", command=self.select_dna_file)
        self.select_file_button.pack(pady=4)

        # Number of clusters entry for KMeans
        self.kmeans_cluster_label = tk.Label(self.root, text="Enter Number of Clusters (KMeans):")
        self.kmeans_cluster_label.pack()
        self.kmeans_cluster_entry = tk.Entry(self.root)
        self.kmeans_cluster_entry.pack(pady=5)

        self.kmeans_button = tk.Button(self.root, text="Run KMeans Clustering", command=self.run_kmeans_clustering)
        self.kmeans_button.pack(pady=5)

        # Number of clusters entry for KMedoids
        self.kmedoids_cluster_label = tk.Label(self.root, text="Enter Number of Clusters (KMedoids):")
        self.kmedoids_cluster_label.pack()
        self.kmedoids_cluster_entry = tk.Entry(self.root)
        self.kmedoids_cluster_entry.pack(pady=5)

        self.kmedoids_button = tk.Button(self.root, text="Run KMedoids Clustering", command=self.run_kmedoids_clustering)
        self.kmedoids_button.pack(pady=5)

        # Number of min. points for DBSCAN
        self.dbscan_min_points_label = tk.Label(self.root, text="Enter Min. Points (DBSCAN):")
        self.dbscan_min_points_label.pack()
        self.dbscan_min_points_entry = tk.Entry(self.root)
        self.dbscan_min_points_entry.pack(pady=5)

        # Neighbourhood radius for DBSCAN
        self.dbscan_radius_label = tk.Label(self.root, text="Enter Neighbourhood Radius (DBSCAN):")
        self.dbscan_radius_label.pack()
        self.dbscan_radius_entry = tk.Entry(self.root)
        self.dbscan_radius_entry.pack(pady=5)

        self.dbscan_button = tk.Button(self.root, text="Run DBSCAN Clustering", command=self.run_dbscan_clustering)
        self.dbscan_button.pack(pady=5)

         # Number of clusters entry for Agglomerative Clustering
        self.agglomerative_cluster_label = tk.Label(self.root, text="Enter Number of Clusters (Agglomerative):")
        self.agglomerative_cluster_label.pack()
        self.agglomerative_cluster_entry = tk.Entry(self.root)
        self.agglomerative_cluster_entry.pack(pady=5)

        self.agglomerative_button = tk.Button(self.root, text="Run Agglomerative Clustering", command=self.run_agglomerative_clustering)
        self.agglomerative_button.pack(pady=5)

        # Show Dendrogram button
        self.show_dendrogram_button = tk.Button(self.root, text="Show Agglomerative Dendrogram", command=self.show_dendrogram)
        self.show_dendrogram_button.pack(pady=5)

        # Result display (Text widget and Scrollbar)
        self.result_text = tk.Text(self.root, wrap=tk.WORD, height=7, width=120)
        self.result_text.pack(pady=5, padx=10)

        self.result_scrollbar = tk.Scrollbar(self.root, command=self.result_text.yview)
        self.result_scrollbar.pack(side=tk.RIGHT, fill=tk.Y)
        self.result_text.config(yscrollcommand=self.result_scrollbar.set)

        # Plot clusters button
        self.plot_button = tk.Button(self.root, text="Plot Clusters", command=self.plot_clusters)
        self.plot_button.pack(pady=5)

        # Silhouette Score label
        self.silhouette_label = tk.Label(self.root, text="Silhouette Score: ")
        self.silhouette_label.pack(pady=5)

    def select_dna_file(self):
        file_path = filedialog.askopenfilename(title="Select DNA File", filetypes=[("Text files", "*.txt")])
        if file_path:
            self.load_dna_sequences(file_path)

    def load_dna_sequences(self, file_path):
        with open(file_path, 'r') as file:
            self.dna_sequences = [line.strip() for line in file]

        # Split data into training and testing sets
        if self.dna_sequences:
            self.dna_sequences_train, self.dna_sequences_test, self.true_labels_train, self.true_labels_test = train_test_split(
                self.dna_sequences, self.true_labels, test_size=0.2, random_state=52
            )

        self.cluster_labels_train = [-1] * len(self.dna_sequences_train)
        self.cluster_labels_test = [-1] * len(self.dna_sequences_test)

    def show_dendrogram(self):
        if not self.dna_sequences:
            self.result_text.insert(tk.END, "No DNA sequences loaded.\n\n")
            return
        # Calculate pairwise Levenshtein distances
        distance_matrix = np.zeros((len(self.dna_sequences), len(self.dna_sequences)))
        for i in range(len(self.dna_sequences)):
            for j in range(i + 1, len(self.dna_sequences)):
                distance_matrix[i, j] = self.calculate_levenshtein_distance(
                    self.dna_sequences[i], self.dna_sequences[j]
                )
                distance_matrix[j, i] = distance_matrix[i, j]

        # Calculate linkage matrix for dendrogram
        linkage_matrix = hierarchy.linkage(distance_matrix, method='average')

        # Create a new Tkinter window for dendrogram
        dendrogram_window = tk.Toplevel(self.root)
        dendrogram_window.title("Dendrogram")

        # Create a figure with adjustable size
        fig, ax = plt.subplots(figsize=(10, 6))

        # Plot dendrogram
        dendrogram = hierarchy.dendrogram(linkage_matrix, ax=ax)

        # Add scrolling capabilities
        canvas = FigureCanvasTkAgg(fig, master=dendrogram_window)
        canvas_widget = canvas.get_tk_widget()
        canvas_widget.pack(side=tk.TOP, fill=tk.BOTH, expand=1)

        # Create a scrollbar for the canvas
        scrollbar = tk.Scrollbar(dendrogram_window, orient=tk.VERTICAL, command=canvas.get_tk_widget().yview)
        scrollbar.pack(side=tk.RIGHT, fill=tk.Y)

        # Configure the canvas to use the scrollbar
        canvas.get_tk_widget().configure(yscrollcommand=scrollbar.set)

        # Embed the dendrogram plot in Tkinter
        canvas.draw()


    def calculate_levenshtein_distance(self, sequence1, sequence2):
        return lev.distance(sequence1, sequence2)

    def run_kmeans_clustering(self):
        num_clusters = self.kmeans_cluster_entry.get()
        try:
            num_clusters = int(num_clusters)
        except ValueError:
            self.result_label.config(text="Invalid input for the number of clusters.")
            return

        Thread(target=self.run_clustering_algorithm, args=(KMeans(n_clusters=num_clusters, random_state=42, n_init=10), "KMeans")).start()

    def run_kmedoids_clustering(self):
        num_clusters = self.kmedoids_cluster_entry.get()

        try:
            num_clusters = int(num_clusters)
        except ValueError:
            self.result_label.config(text="Invalid input for the number of clusters.")
            return

        # Use threading to run the clustering algorithm
        Thread(target=self.run_clustering_algorithm, args=(KMedoids(n_clusters=num_clusters, random_state=42), "KMedoids")).start()

    def run_dbscan_clustering(self):
        min_points = self.dbscan_min_points_entry.get()
        radius = self.dbscan_radius_entry.get()

        try:
            min_points = int(min_points)
            radius = float(radius)
        except ValueError:
            self.result_label.config(text="Invalid input for DBSCAN parameters.")
            return

        Thread(target=self.run_clustering_algorithm, args=(DBSCAN(min_samples=min_points, eps=radius), "DBSCAN")).start()

    def run_agglomerative_clustering(self):
        num_clusters = self.agglomerative_cluster_entry.get()

        try:
            num_clusters = int(num_clusters)
        except ValueError:
            self.result_label.config(text="Invalid input for the number of clusters.")
            return

        Thread(target=self.run_clustering_algorithm, args=(AgglomerativeClustering(n_clusters=num_clusters), "Agglomerative Hierarchical")).start()


    def run_clustering_algorithm(self, algorithm, algorithm_name):
        if not self.dna_sequences:
            self.result_label.config(text="No DNA sequences loaded.")
            return

        # Calculate pairwise Levenshtein distances
        distance_matrix = np.zeros((len(self.dna_sequences), len(self.dna_sequences)))
        for i in range(len(self.dna_sequences)):
            for j in range(i + 1, len(self.dna_sequences)):
                distance_matrix[i, j] = self.calculate_levenshtein_distance(
                    self.dna_sequences[i], self.dna_sequences[j]
                )
                distance_matrix[j, i] = distance_matrix[i, j]

        # Apply clustering algorithm
        clusters = algorithm.fit_predict(distance_matrix)

        self.cluster_labels = clusters

        # Check if all points are noise (cluster label is -1 for DBSCAN)
        if algorithm_name == "DBSCAN" and np.all(clusters == -1):
            self.result_label.config(text="DBSCAN identified all points as noise.")
            return

        # Calculate silhouette score
        silhouette_avg = silhouette_score(distance_matrix, clusters)

        # Display clustering results in the Text widget
        result_str = f"{algorithm_name} Clustering Results:\nSilhouette Score: {silhouette_avg:.2f}\n"
        for cluster_id in np.unique(clusters):
            cluster_indices = np.where(clusters == cluster_id)[0]
            if cluster_id == -1 and algorithm_name == "DBSCAN":
                result_str += f"Noise Points: {', '.join(map(str, cluster_indices))}\n"
            else:
                result_str += f"Cluster {cluster_id + 1}: {', '.join(map(str, cluster_indices))}\n"

        self.result_text.insert(tk.END, result_str + "\n")
        # Update the Silhouette Score label in the GUI
        self.silhouette_label.config(text=f"Silhouette Score: {silhouette_avg:.2f}")


    def plot_clusters(self):
        if not self.dna_sequences or len(self.cluster_labels) == 0:
            return

        unique_clusters = np.unique(self.cluster_labels)

        # Create a new Tkinter window
        plot_window = tk.Toplevel(self.root)
        plot_window.title("Cluster Plot")

        # Create a figure and axis with a specified size
        fig, ax = plt.subplots(figsize=(15, 8))

        # Plot points for each cluster
        for cluster_id in unique_clusters:
            if cluster_id == -1:
                # Plot noise points separately
                noise_indices = np.where(self.cluster_labels == -1)[0]
                ax.scatter([-1] * len(noise_indices), noise_indices, label='Noise', marker='x')
            else:
                cluster_indices = np.where(self.cluster_labels == cluster_id)[0]
                cluster_dna_numbers = [i for i in cluster_indices]  # Keep it 0-based index
                ax.scatter([cluster_id] * len(cluster_indices), cluster_dna_numbers, label=f'Cluster {cluster_id + 1}')

        # Set y-axis ticks to whole numbers starting from 0 and display ordinal labels
        ax.set_yticks(range(len(self.dna_sequences)))

        p = inflect.engine()
        ax.set_yticklabels([p.ordinal(i) for i in range(len(self.dna_sequences))])

        ax.set_xlabel('Cluster')
        ax.set_ylabel('DNA Sequences')
        ax.legend()

        # Remove previous plot if it exists
        for widget in self.root.winfo_children():
            if isinstance(widget, FigureCanvasTkAgg):
                widget.get_tk_widget().destroy()

        # Embed the plot in Tkinter
        canvas = FigureCanvasTkAgg(fig, master=plot_window)
        canvas_widget = canvas.get_tk_widget()
        canvas_widget.pack()

        canvas.draw()

if __name__ == "__main__":
    root = tk.Tk()
    app = DNAClusteringApp(root)
    root.mainloop()
