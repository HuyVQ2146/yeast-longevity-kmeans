import pandas as pd
from math import sqrt
import random

def normalization(value):
    current_max = 10.0
    current_min = -10.0

    new_max = 1.0
    new_min = 0
    
    return ((value - current_min) / (current_max - current_min)) * (new_max - new_min) + new_min

def euclid_distance(num1, num2, num3):
    return sqrt(num1**2 + num2**2 + num3**2)

def k_means_clustering(k, max_iterations, data, centroids):
    for _ in range(max_iterations):
        clusters = {i: [] for i in range(k)}

        for idx, row in data.iterrows():
            name = row["Name"]
            x = row["euclid_distance"]
            known = row["Known"]

            distances = [abs(x - c) for c in centroids]
            nearest = distances.index(min(distances))

            clusters[nearest].append((name, x, known))

        for i in range(k):
            if len(clusters[i]) > 0:
                centroids[i] = sum(item[1] for item in clusters[i]) / len(clusters[i])
    return clusters

def calculate_silhouette_score(clusters):
    all_scores = []

    for cluster_idx, points in clusters.items():
        n_points = len(points)
        if n_points == 0:
            continue

        for name, value, known in points:
            if n_points == 1:
                a_i = 0.0
            else:
                dist_to_others = sum(
                    abs(value - other_val)
                    for other_name, other_val, _ in points
                    if other_name != name
                )
                a_i = dist_to_others / (n_points - 1)

            nearest_other_avg = float("inf")
            has_other_cluster = False

            for other_idx, other_points in clusters.items():
                if other_idx == cluster_idx or len(other_points) == 0:
                    continue
                has_other_cluster = True
                avg_dist = sum(
                    abs(value - other_val)
                    for _, other_val, _ in other_points
                ) / len(other_points)
                nearest_other_avg = min(nearest_other_avg, avg_dist)

            b_i = nearest_other_avg if has_other_cluster else 0.0

            denom = max(a_i, b_i)
            s_i = (b_i - a_i) / denom if denom != 0 else 0.0

            all_scores.append(s_i)

    return sum(all_scores) / len(all_scores) if all_scores else 0.0


def main():
    data = pd.read_csv('cleaned_data.csv')
    
    data["normalized_1"] = data["sch9/wt"].apply(normalization)
    data["normalized_2"] = data["ras2/wt"].apply(normalization)
    data["normalized_3"] = data["tor1/wt"].apply(normalization)
    data["euclid_distance"] = data.apply(lambda row: euclid_distance(row["normalized_1"], row["normalized_2"], row["normalized_3"]), axis=1)
    
    random.seed(42)
    X = data["euclid_distance"].tolist()

    print("=" * 50)
    print("Silhouette Score for K = 2 to 10")
    print("=" * 50)

    results = []

    for test_k in range(2, 11):
        random.seed(42)
        initial_centroids = random.sample(X, test_k)

        clusters = k_means_clustering(test_k, 100, data, initial_centroids)

        score = calculate_silhouette_score(clusters)

        cluster_sizes = [len(clusters[i]) for i in range(test_k)]
        results.append((test_k, score, cluster_sizes))

        print(f"K = {test_k:2d}  |  Silhouette Score = {score:.4f}  |  Cluster sizes = {cluster_sizes}")

    best_k, best_score, _ = max(results, key=lambda r: r[1])
    print("-" * 50)
    print(f"Best K = {best_k} with Silhouette Score = {best_score:.4f}")
    
    random.seed(42)
    initial_centroids = random.sample(X, best_k)
    best_clusters = k_means_clustering(best_k, 100, data, initial_centroids)
    score = calculate_silhouette_score(best_clusters)
    cluster_sizes = [len(best_clusters[i]) for i in range(best_k)]
    results.append((best_k, score, cluster_sizes))

    print("=" * 65)
    print(f"{'Cluster':<10} {'Total':<8} {'Known=1':<8} {'Known=0':<8} {'% Known':<10}")
    print("-" * 65)

    total_known = 0
    total_genes = 0

    for cluster_idx in range(best_k):
        points = best_clusters[cluster_idx]
        n_total = len(points)
        n_known = sum(1 for _, _, known in points if known == 1)
        n_unknown = n_total - n_known
        pct_known = (n_known / n_total * 100) if n_total > 0 else 0.0
        total_known += n_known
        total_genes += n_total
        print(f"{cluster_idx+1:<10} {n_total:<8} {n_known:<8} {n_unknown:<8} {pct_known:>6.1f}%")

    print("-" * 65)
    overall_pct = (total_known / total_genes * 100) if total_genes > 0 else 0.0
    print(f"{'TOTAL':<10} {total_genes:<8} {total_known:<8} {total_genes - total_known:<8} {overall_pct:>6.1f}%")


    print()
    print("=" * 65)
    print("Gene Names in Each Cluster")
    print("=" * 65)

    for cluster_idx in range(best_k):
        points = best_clusters[cluster_idx]
        names = [name for name, _, _ in points]
        known_names = [name for name, _, known in points if known == 1]
        candidate_names = [name for name, _, known in points if known == 0]

        print(f"\n--- Cluster {cluster_idx+1} ({len(points)} genes) ---")
        print(f"  All genes     : {', '.join(names)}")
        if known_names:
            print(f"  Known genes   : {', '.join(known_names)}")
        else:
            print(f"  Known genes   : (none)")
        if candidate_names:
            print(f"  Candidate genes: {', '.join(candidate_names)}")
        else:
            print(f"  Candidate genes: (none)")


    
if __name__ == "__main__":
    main()
