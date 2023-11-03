import numpy as np
from IPython.display import display, Math

def create_measurement_covariance_matrix(diagonal_terms):
    """
    Creates a square covariance matrix known as Sigma_{I_{data}}.
    
    Parameters:
    diagonal_terms (list): The values for the diagonal of the matrix.
    
    Returns:
    numpy.ndarray: The covariance matrix with the given diagonal terms.
    """
    n = len(diagonal_terms)
    matrix = np.zeros((n, n))
    np.fill_diagonal(matrix, diagonal_terms)
    return matrix

def construct_QU_data(num_standards, Qs, Us):
    """
    Constructs the QU_data matrix with entire columns filled with Q's and U's.
    
    Parameters:
    num_standards (int): The number of standards.
    Qs (list): The list of Q values.
    Us (list): The list of U values.
    
    Returns:
    numpy.ndarray: The QU_data matrix.
    """
    if len(Qs) != 2 * num_standards or len(Us) != 2 * num_standards:
        raise ValueError("The length of Q's and U's must match two times the number of standards.")
    
    QU_data = np.zeros((num_standards * 2, 2))
    QU_data[:, 0] = Qs  # Fill the entire first column with Q's
    QU_data[:, 1] = Us  # Fill the entire second column with U's
    
    return QU_data


def calculate_instrument_param_errs(Sigma_I_data, QU_data):
    """
    Calculates the covariance matrix on the instrument parameters, referred to as
    Sigma_{eta, chi}.
    
    Parameters:
    Sigma_I_data (numpy.ndarray): The measurement covariance matrix.
    QU_data (numpy.ndarray): The QU_data matrix.
    
    Returns:
    numpy.ndarray: The covariance matrix of instrument parameters.
    """
    QU_data_transpose = QU_data.T
    Sigma_I_data_inv = np.linalg.inv(Sigma_I_data)
    intermediate = QU_data_transpose @ Sigma_I_data_inv @ QU_data
    Sigma_eta_xi = np.linalg.inv(intermediate)
    
    # Display the diagonal terms of the Sigma_eta_xi matrix in LaTeX format
    display(Math(f'\\text{{Error on }} \\eta_{{Q}}: {Sigma_eta_xi[0, 0]:.2f}'))
    display(Math(f'\\text{{Error on }} \\chi_{{U \\to Q}}: {Sigma_eta_xi[1, 1]:.2f}'))
    
    return Sigma_eta_xi
