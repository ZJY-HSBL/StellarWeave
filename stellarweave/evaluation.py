"""One-to-one positional matching and synthetic detection metrics."""
import numpy as np
from scipy.optimize import linear_sum_assignment
from scipy.spatial.distance import cdist

def match_sources(truth_xy, predicted_xy, radius=2.0):
    truth = np.asarray(truth_xy, dtype=float).reshape(-1, 2)
    pred = np.asarray(predicted_xy, dtype=float).reshape(-1, 2)
    if not len(truth) or not len(pred):
        return np.empty((0, 3), dtype=float)
    d = cdist(truth, pred)
    n, m = d.shape
    penalty = float(radius) + .01
    cost = np.full((n + m, n + m), 1e5)
    cost[:n, :m] = np.where(d <= radius, d, 1e5)
    cost[np.arange(n), m + np.arange(n)] = penalty
    cost[n + np.arange(m), np.arange(m)] = penalty
    cost[n:, m:] = 0
    rows, cols = linear_sum_assignment(cost)
    pairs = [(i, j, d[i, j]) for i, j in zip(rows, cols)
             if i < n and j < m and d[i, j] <= radius]
    return np.asarray(pairs, dtype=float).reshape(-1, 3)

def detection_metrics(truth_xy, predicted_xy, radius=2.0):
    n, m = len(truth_xy), len(predicted_xy)
    pairs = match_sources(truth_xy, predicted_xy, radius)
    tp = len(pairs)
    return {"true_sources": n, "detected_sources": m, "matched": tp,
            "false_positive": m - tp, "missed": n - tp,
            "precision": tp / m if m else 0.0,
            "recall": tp / n if n else 0.0,
            "f1": 2 * tp / (n + m) if n + m else 0.0,
            "mean_pos_error_px": float(pairs[:, 2].mean()) if tp else None}

def flux_binned_recall(truth, predicted_xy,
                       flux_bins=(0, 300, 800, 2000, 6000, np.inf), radius=2.0):
    truth = np.asarray(truth)
    pairs = match_sources(truth[:, :2], predicted_xy, radius)
    matched = set(pairs[:, 0].astype(int)) if len(pairs) else set()
    rows = []
    for low, high in zip(flux_bins[:-1], flux_bins[1:]):
        indices = np.flatnonzero((truth[:, 2] >= low) & (truth[:, 2] < high))
        count = sum(int(i) in matched for i in indices)
        rows.append({"flux_min": float(low),
                     "flux_max": float(high) if np.isfinite(high) else None,
                     "total": len(indices), "recovered": count,
                     "recall": count / len(indices) if len(indices) else None})
    return rows
