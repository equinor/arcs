from arcs.traversal import (
    filter_out_compounds_under_probability_threshold,
    calculate_compound_probabilities_relative_to_compound_values,
)


def test_probability_is_evenly_distributed_across_equal_concs():
    concs = {"a": 1, "b": 1}

    probabilities = calculate_compound_probabilities_relative_to_compound_values(concs)

    assert probabilities == {"a": 0.5, "b": 0.5}


def test_probabilities_over_threshold_is_filtered_out():
    probabilities = {"a": 0.1, "b": 0.8, "c": 0.1}
    probability_threshold = 0.1

    filtered_probabilities = (
        filter_out_compounds_under_probability_threshold(
            probability_threshold=probability_threshold, probabilities=probabilities
        )
    )

    assert filtered_probabilities == {"b": 0.8}


def test_probabilities_close_to_threshold_is_kept():
    probability_threshold = 0.1
    probabilities = {"a": probability_threshold + 1e-10, "b": 0.9}

    filtered_probabilities = (
        filter_out_compounds_under_probability_threshold(
            probability_threshold=probability_threshold, probabilities=probabilities
        )
    )

    assert filtered_probabilities == probabilities
