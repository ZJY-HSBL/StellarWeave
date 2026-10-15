import numpy as np
from stellarweave.simulation import simulated_field
from stellarweave.detection import robust_background, scipy_peak_detector
from stellarweave.evaluation import match_sources, detection_metrics, flux_binned_recall

def test_reproducibility():
    image_a, truth_a = simulated_field(size=64, n_stars=5, seed=2)
    image_b, truth_b = simulated_field(size=64, n_stars=5, seed=2)
    np.testing.assert_array_equal(image_a, image_b)
    np.testing.assert_array_equal(truth_a, truth_b)

def test_perfect_match():
    points = np.array([[4., 4.], [9., 9.]])
    assert detection_metrics(points, points)["f1"] == 1
    assert len(match_sources(points, points)) == 2

def test_unmatched_false_positive():
    result = detection_metrics(np.array([[4., 4.]]),
                               np.array([[4., 4.], [30., 30.]]))
    assert result["false_positive"] == 1
    assert result["precision"] == .5

def test_one_to_one():
    assert len(match_sources([[5., 5.], [5.2, 5.]], [[5.1, 5.]])) == 1

def test_synthetic_detection():
    image, truth = simulated_field(size=160, n_stars=15, seed=7, min_flux=1000)
    found = scipy_peak_detector(image, nsigma=4)
    assert len(found) > 5
    assert detection_metrics(truth[:, :2], found)["recall"] > .3

def test_background():
    image = np.full((20, 20), 13.)
    image[4, 4] = 300
    median, _ = robust_background(image)
    assert abs(median - 13) < 1

def test_flux_binning():
    truth = np.array([[4., 4., 100.], [20., 20., 900.]])
    bins = flux_binned_recall(truth, [[4., 4.]], flux_bins=(0, 500, 1000))
    assert bins[0]["recall"] == 1
    assert bins[1]["recall"] == 0
