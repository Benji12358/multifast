module boundaries_types
    implicit none

    enum, bind(c)
    enumerator :: periodic, Dirichlet, symetric, antisymetric, values, values_and_derivatives

end enum


end module boundaries_types
