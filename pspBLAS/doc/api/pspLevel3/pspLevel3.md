# Overview

The `pspLevel3.F90` file defines the `pspLevel3` module, which acts as a central aggregator for various Level 3 Basic Linear Algebra Subprograms (BLAS) operations within the pspBLAS library. Level 3 BLAS operations are characterized by matrix-matrix computations, which typically offer the highest potential for data reuse and performance on modern computer architectures.

This module does not directly implement computational routines itself. Instead, it serves to group and re-export functionalities from other, more specialized modules that contain the actual implementations of different types of matrix-matrix operations.

# Key Components

*   **Module `pspLevel3`**:
    *   This is the primary entity defined in the file. Its role is to consolidate various Level 3 operation modules from the pspBLAS library. By including a `USE pspLevel3` statement (often indirectly via `USE pspBLAS`), users gain access to the collective public interfaces of these specialized modules.

# Important Variables/Constants

*   `HAVE_CONFIG_H`: A preprocessor macro. If defined, it includes a `config.h` file. This file is typically used in C/C++ and sometimes Fortran projects to manage compile-time configurations and preprocessor definitions that might affect how the library is built or behaves (e.g., enabling or disabling certain features, setting paths, etc.).

# Usage Examples

Users of the pspBLAS library would typically not `USE pspLevel3` directly but would rather `USE pspBLAS`, which in turn makes the functionalities aggregated by `pspLevel3` available.

```fortran
USE pspBLAS  ! This statement also provides access to pspLevel3 contents

! After this, users can call various Level 3 routines, such as:
! - General matrix-matrix multiplication (from pspGemm module)
!   CALL psp_gemm(M, N, K, A_local, opA, B_local, opB, C_local, alpha, beta)
!
! - Sparse matrix-dense matrix multiplication (conceptually, from pspSpmm)
!   CALL psp_spmm(...)
!
! - Dense matrix-sparse matrix multiplication (conceptually, from pspMspm)
!   CALL psp_mspm(...)
!
! - Sparse matrix-sparse matrix multiplication (conceptually, from pspSpmSpm)
!   CALL psp_spspmm(...)
!
! - Matrix summation/addition (conceptually, from pspMatSum)
!   CALL psp_matsum(...)

! Note: The exact names and signatures of the routines depend on their definitions
! in the respective specialized modules.
```

# Dependencies and Interactions

The `pspLevel3` module is an aggregator and thus has dependencies on several other modules that provide specific Level 3 functionalities:

*   **`pspMatSum`**: This module is `USE`d, suggesting it provides routines for matrix addition or summation operations (e.g., `C = alpha*A + beta*B`).
*   **`pspGemm`**: This module is `USE`d and contains implementations for general matrix-matrix multiplication (GEMM), typically for dense matrices, such as `C = alpha*op(A)*op(B) + beta*C`.
*   **`pspSpmm`**: This module is `USE`d, likely providing routines for Sparse Matrix - Dense Matrix Multiplication (SpMM), where one of the input matrices is sparse and the other is dense.
*   **`pspMspm`**: This module is `USE`d. The naming suggests it might handle Dense Matrix - Sparse Matrix Multiplication (MSpM), which is functionally similar to SpMM but might have different algorithmic approaches or interfaces.
*   **`pspSpmSpm`**: This module is `USE`d, indicating it provides routines for Sparse Matrix - Sparse Matrix Multiplication (SpSpMM), where both input matrices are in a sparse format.

The `pspLevel3` module itself is a component of the larger `pspBLAS` library, as the main `pspBLAS` module includes `USE pspLevel3`. This hierarchical structure helps organize the library and allows users to access a wide range of matrix operations through a single top-level module.
