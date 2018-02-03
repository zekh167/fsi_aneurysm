!!> \file
!> $Id: Stokes_equations_routines.f90 372 2009-04-20
!> \author Sebastian Krittian
!> \brief This module handles all Stokes fluid routines.
!>
!> \section LICENSE
!>
!> Version: MPL 1.1/GPL 2.0/LGPL 2.1
!>
!> The contents of this file are subject to the Mozilla Public License
!> Version 1.1 (the "License"); you may not use this file except in
!> compliance with the License. You may obtain a copy of the License at
!> http://www.mozilla.org/MPL/
!>
!> Software distributed under the License is distributed on an "AS IS"
!> basis, WITHOUT WARRANTY OF ANY KIND, either express or implied. See the
!> License for the specific language governing rights and limitations
!> under the License.
!>
!> The Original Code is OpenCMISS
!>
!> The Initial Developer of the Original Code is University of Auckland,
!> Auckland, New Zealand, the University of Oxford, Oxford, United
!> Kingdom and King's College, London, United Kingdom. Portions created
!> by the University of Auckland, the University of Oxford and King's
!> College, London are Copyright (C) 2007-2010 by the University of
!> Auckland, the University of Oxford and King's College, London.
!> All Rights Reserved.
!>
!> Contributor(s): Sebastian Krittian, Chris Bradley
!>
!> Alternatively, the contents of this file may be used under the terms of
!> either the GNU General Public License Version 2 or later (the "GPL"), or
!> the GNU Lesser General Public License Version 2.1 or later (the "LGPL"),
!> in which case the provisions of the GPL or the LGPL are applicable instead
!> of those above. If you wish to allow use of your version of this file only
!> under the terms of either the GPL or the LGPL, and not to allow others to
!> use your version of this file under the terms of the MPL, indicate your
!> decision by deleting the provisions above and replace them with the notice
!> and other provisions required by the GPL or the LGPL. If you do not delete
!> the provisions above, a recipient may use your version of this file under
!> the terms of any one of the MPL, the GPL or the LGPL.
!>

!>This module handles all Stokes fluid routines.
MODULE STOKES_EQUATIONS_ROUTINES

  USE ANALYTIC_ANALYSIS_ROUTINES
  USE BaseRoutines
  USE BASIS_ROUTINES
  USE BOUNDARY_CONDITIONS_ROUTINES
  USE Constants
  USE CONTROL_LOOP_ROUTINES
  USE ControlLoopAccessRoutines
  USE DISTRIBUTED_MATRIX_VECTOR
  USE DOMAIN_MAPPINGS
  USE EquationsRoutines
  USE EquationsAccessRoutines
  USE EquationsMappingRoutines
  USE EquationsMatricesRoutines
  USE EQUATIONS_SET_CONSTANTS
  USE EquationsSetAccessRoutines
  USE FIELD_ROUTINES
  USE FIELD_IO_ROUTINES
  USE FieldAccessRoutines
  USE FLUID_MECHANICS_IO_ROUTINES
  USE INPUT_OUTPUT
  USE ISO_VARYING_STRING
  USE Kinds
  USE MATRIX_VECTOR
  USE NODE_ROUTINES
  USE PROBLEM_CONSTANTS
  USE Strings
  USE SOLVER_ROUTINES
  USE SolverAccessRoutines
  USE Timer
  USE Types

#include "macros.h"

  IMPLICIT NONE

  PRIVATE

  PUBLIC  Stokes_EquationsSetSpecificationSet

  PUBLIC  Stokes_EquationsSetSolutionMethodSet

  PUBLIC  STOKES_EQUATIONS_SET_SETUP

  PUBLIC  Stokes_BoundaryConditionsAnalyticCalculate

  PUBLIC  Stokes_ProblemSpecificationSet

  PUBLIC  STOKES_PROBLEM_SETUP

  PUBLIC  STOKES_FINITE_ELEMENT_CALCULATE

  PUBLIC  STOKES_POST_SOLVE

  PUBLIC  STOKES_PRE_SOLVE

  PUBLIC  STOKES_EQUATION_ANALYTIC_FUNCTIONS

CONTAINS

!
!================================================================================================================================
!

  !>Sets/changes the solution method for a Stokes flow equation type of an fluid mechanics equations set class.
  SUBROUTINE Stokes_EquationsSetSolutionMethodSet(EQUATIONS_SET,SOLUTION_METHOD,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to set the solution method for
    INTEGER(INTG), INTENT(IN) :: SOLUTION_METHOD !<The solution method to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("Stokes_EquationsSetSolutionMethodSet",err,error,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(.NOT.ALLOCATED(EQUATIONS_SET%SPECIFICATION)) THEN
        CALL FlagError("Equations set specification is not allocated.",err,error,*999)
      ELSE IF(SIZE(EQUATIONS_SET%SPECIFICATION,1)/=3) THEN
        CALL FlagError("Equations set specification must have 3 entries for a Stokes equations set.", &
          & err,error,*999)
      END IF
      SELECT CASE(EQUATIONS_SET%SPECIFICATION(3))
      CASE(EQUATIONS_SET_STATIC_STOKES_SUBTYPE, &
        & EQUATIONS_SET_LAPLACE_STOKES_SUBTYPE, &
        & EQUATIONS_SET_TRANSIENT_STOKES_SUBTYPE, &
        & EQUATIONS_SET_ALE_STOKES_SUBTYPE, &
        & EQUATIONS_SET_PGM_STOKES_SUBTYPE)
        SELECT CASE(SOLUTION_METHOD)
          CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
            EQUATIONS_SET%SOLUTION_METHOD=EQUATIONS_SET_FEM_SOLUTION_METHOD
          CASE(EQUATIONS_SET_BEM_SOLUTION_METHOD)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE(EQUATIONS_SET_FD_SOLUTION_METHOD)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE(EQUATIONS_SET_FV_SOLUTION_METHOD)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE(EQUATIONS_SET_GFEM_SOLUTION_METHOD)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE(EQUATIONS_SET_GFV_SOLUTION_METHOD)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE DEFAULT
            localError="The specified solution method of "//TRIM(NumberToVString(SOLUTION_METHOD,"*",err,error))//" is invalid."
            CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE DEFAULT
        localError="The third equations set specification of "// &
          & TRIM(NumberToVstring(EQUATIONS_SET%SPECIFICATION(3),"*",err,error))// &
          & " is not valid for a Stokes flow equation of a fluid mechanics equations set."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    ELSE
      CALL FlagError("Equations set is not associated.",err,error,*999)
    ENDIF

    EXITS("Stokes_EquationsSetSolutionMethodSet")
    RETURN
999 ERRORSEXITS("Stokes_EquationsSetSolutionMethodSet",err,error)
    RETURN 1

  END SUBROUTINE Stokes_EquationsSetSolutionMethodSet

!
!================================================================================================================================
!

  !>Sets the equation specification for a Stokes flow equation of a fluid mechanics equations set.
  SUBROUTINE Stokes_EquationsSetSpecificationSet(equationsSet,specification,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet !<A pointer to the equations set to set the specification for
    INTEGER(INTG), INTENT(IN) :: specification(:) !<The equations set specification to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
    INTEGER(INTG) :: subtype

    ENTERS("Stokes_EquationsSetSpecificationSet",err,error,*999)

    IF(ASSOCIATED(equationsSet)) THEN
      IF(SIZE(specification,1)/=3) THEN
        CALL FlagError("Equations set specification must have three entries for a Stokes type equations set.", &
          & err,error,*999)
      END IF
      subtype=specification(3)
      SELECT CASE(subtype)
      CASE(EQUATIONS_SET_STATIC_STOKES_SUBTYPE, &
          & EQUATIONS_SET_LAPLACE_STOKES_SUBTYPE, &
          & EQUATIONS_SET_TRANSIENT_STOKES_SUBTYPE, &
          & EQUATIONS_SET_ALE_STOKES_SUBTYPE, &
          & EQUATIONS_SET_PGM_STOKES_SUBTYPE)
        !ok
      CASE(EQUATIONS_SET_OPTIMISED_STOKES_SUBTYPE)
        CALL FlagError("Not implemented yet.",err,error,*999)
      CASE DEFAULT
        localError="The third equations set specification of "//TRIM(NumberToVstring(specification(3),"*",err,error))// &
          & " is not valid for Stokes flow of a fluid mechanics equations set."
        CALL FlagError(localError,err,error,*999)
      END SELECT
      !Set full specification
      IF(ALLOCATED(equationsSet%specification)) THEN
        CALL FlagError("Equations set specification is already allocated.",err,error,*999)
      ELSE
        ALLOCATE(equationsSet%specification(3),stat=err)
        IF(err/=0) CALL FlagError("Could not allocate equations set specification.",err,error,*999)
      END IF
      equationsSet%specification(1:3)=[EQUATIONS_SET_FLUID_MECHANICS_CLASS,EQUATIONS_SET_STOKES_EQUATION_TYPE,subtype]
    ELSE
      CALL FlagError("Equations set is not associated.",err,error,*999)
    END IF

    EXITS("Stokes_EquationsSetSpecificationSet")
    RETURN
999 ERRORS("Stokes_EquationsSetSpecificationSet",err,error)
    EXITS("Stokes_EquationsSetSpecificationSet")
    RETURN 1

  END SUBROUTINE Stokes_EquationsSetSpecificationSet

!
!================================================================================================================================
!

  !>Sets up the standard Stokes fluid setup.
  SUBROUTINE STOKES_EQUATIONS_SET_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to setup
    TYPE(EQUATIONS_SET_SETUP_TYPE), INTENT(INOUT) :: EQUATIONS_SET_SETUP !<The equations set setup information
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: GEOMETRIC_SCALING_TYPE,GEOMETRIC_MESH_COMPONENT
    TYPE(DECOMPOSITION_TYPE), POINTER :: GEOMETRIC_DECOMPOSITION
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(EQUATIONS_SET_MATERIALS_TYPE), POINTER :: EQUATIONS_MATERIALS
    TYPE(VARYING_STRING) :: localError
    INTEGER(INTG):: DEPENDENT_FIELD_NUMBER_OF_VARIABLES,DEPENDENT_FIELD_NUMBER_OF_COMPONENTS
    INTEGER(INTG):: INDEPENDENT_FIELD_NUMBER_OF_VARIABLES,INDEPENDENT_FIELD_NUMBER_OF_COMPONENTS
    INTEGER(INTG):: NUMBER_OF_DIMENSIONS,GEOMETRIC_COMPONENT_NUMBER
    INTEGER(INTG):: MATERIAL_FIELD_NUMBER_OF_VARIABLES,MATERIAL_FIELD_NUMBER_OF_COMPONENTS,I

    ENTERS("STOKES_EQUATIONS_SET_SETUP",err,error,*999)

    NULLIFY(EQUATIONS)
    NULLIFY(vectorMapping)
    NULLIFY(vectorMatrices)
    NULLIFY(GEOMETRIC_DECOMPOSITION)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(.NOT.ALLOCATED(EQUATIONS_SET%SPECIFICATION)) THEN
        CALL FlagError("Equations set specification is not allocated.",err,error,*999)
      ELSE IF(SIZE(EQUATIONS_SET%SPECIFICATION,1)<3) THEN
        CALL FlagError("Equations set specification must have >= 3 entries for a Stokes flow equations set.", &
          & err,error,*999)
      END IF
      SELECT CASE(EQUATIONS_SET%SPECIFICATION(3))
        !Select Stokes subtypes
        CASE(EQUATIONS_SET_STATIC_STOKES_SUBTYPE, &
          & EQUATIONS_SET_TRANSIENT_STOKES_SUBTYPE, &
          & EQUATIONS_SET_LAPLACE_STOKES_SUBTYPE, &
          & EQUATIONS_SET_ALE_STOKES_SUBTYPE, &
          & EQUATIONS_SET_PGM_STOKES_SUBTYPE)
          SELECT CASE(EQUATIONS_SET_SETUP%SETUP_TYPE)
            !Set solution method
            CASE(EQUATIONS_SET_SETUP_INITIAL_TYPE)
              SELECT CASE(EQUATIONS_SET%SPECIFICATION(3))
                CASE(EQUATIONS_SET_STATIC_STOKES_SUBTYPE,EQUATIONS_SET_LAPLACE_STOKES_SUBTYPE, &
                  & EQUATIONS_SET_TRANSIENT_STOKES_SUBTYPE,EQUATIONS_SET_ALE_STOKES_SUBTYPE, &
                  & EQUATIONS_SET_PGM_STOKES_SUBTYPE)
                  SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
                    CASE(EQUATIONS_SET_SETUP_START_ACTION)
                      CALL Stokes_EquationsSetSolutionMethodSet(EQUATIONS_SET ,&
                        & EQUATIONS_SET_FEM_SOLUTION_METHOD,err,error,*999)
                      EQUATIONS_SET%SOLUTION_METHOD=EQUATIONS_SET_FEM_SOLUTION_METHOD
                      CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
                      !!TODO: Check valid setup
                    CASE DEFAULT
                      localError="The action type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%ACTION_TYPE, &
                        & "*",err,error))// " for a setup type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP% &
                        & SETUP_TYPE,"*",err,error))// " is invalid for a standard Stokes fluid."
                      CALL FlagError(localError,err,error,*999)
                  END SELECT
                CASE DEFAULT
                  localError="The third equations set specification of "// &
                    & TRIM(NumberToVstring(EQUATIONS_SET%SPECIFICATION(3),"*", &
                    & err,error))//" is invalid for a Stokes flow equations set."
                  CALL FlagError(localError,err,error,*999)
              END SELECT
            CASE(EQUATIONS_SET_SETUP_GEOMETRY_TYPE)
              !Set geometric field
              SELECT CASE(EQUATIONS_SET%SPECIFICATION(3))
                CASE(EQUATIONS_SET_STATIC_STOKES_SUBTYPE,EQUATIONS_SET_LAPLACE_STOKES_SUBTYPE, &
                  & EQUATIONS_SET_TRANSIENT_STOKES_SUBTYPE,EQUATIONS_SET_ALE_STOKES_SUBTYPE, &
                  & EQUATIONS_SET_PGM_STOKES_SUBTYPE)
                    !Do nothing???
                CASE DEFAULT
                  localError="The third equations set specification of "// &
                    & TRIM(NumberToVString(EQUATIONS_SET%SPECIFICATION(3),"*", &
                    & err,error))//" is invalid for a Stokes flow equations set."
                  CALL FlagError(localError,err,error,*999)
              END SELECT
            CASE(EQUATIONS_SET_SETUP_DEPENDENT_TYPE)
              !Set dependent field
              SELECT CASE(EQUATIONS_SET%SPECIFICATION(3))
                CASE(EQUATIONS_SET_STATIC_STOKES_SUBTYPE,EQUATIONS_SET_LAPLACE_STOKES_SUBTYPE, &
                  & EQUATIONS_SET_TRANSIENT_STOKES_SUBTYPE,EQUATIONS_SET_ALE_STOKES_SUBTYPE, &
                  & EQUATIONS_SET_PGM_STOKES_SUBTYPE)
                  SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
                   !Set start action
                    CASE(EQUATIONS_SET_SETUP_START_ACTION)
                      IF(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD_AUTO_CREATED) THEN
                        !Create the auto created dependent field
                        !start field creation with name 'DEPENDENT_FIELD'
                        CALL FIELD_CREATE_START(EQUATIONS_SET_SETUP%FIELD_USER_NUMBER,EQUATIONS_SET%REGION, &
                          & EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,err,error,*999)
                        !start creation of a new field
                        CALL FIELD_TYPE_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_GENERAL_TYPE,err,error,*999)
                        !label the field
                        CALL FIELD_LABEL_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,"U",err,error,*999)
                        !define new created field to be dependent
                        CALL FIELD_DEPENDENT_TYPE_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                          & FIELD_DEPENDENT_TYPE,err,error,*999)
                        !look for decomposition rule already defined
                        CALL FIELD_MESH_DECOMPOSITION_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_DECOMPOSITION, &
                          & err,error,*999)
                        !apply decomposition rule found on new created field
                        CALL FIELD_MESH_DECOMPOSITION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                          & GEOMETRIC_DECOMPOSITION,err,error,*999)
                        !point new field to geometric field
                        CALL FIELD_GEOMETRIC_FIELD_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,EQUATIONS_SET%GEOMETRY% &
                          & GEOMETRIC_FIELD,err,error,*999)
                        !set number of variables to 2 (1 for U and one for DELUDELN)
                        DEPENDENT_FIELD_NUMBER_OF_VARIABLES=2
                        CALL FIELD_NUMBER_OF_VARIABLES_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                          & DEPENDENT_FIELD_NUMBER_OF_VARIABLES,err,error,*999)
                        CALL FIELD_VARIABLE_TYPES_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,[FIELD_U_VARIABLE_TYPE, &
                          & FIELD_DELUDELN_VARIABLE_TYPE],err,error,*999)
                        CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                          & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
                        CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                          & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
                        CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                          & FIELD_DP_TYPE,err,error,*999)
                        CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                          & FIELD_DP_TYPE,err,error,*999)
                        CALL FIELD_NUMBER_OF_COMPONENTS_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                          & NUMBER_OF_DIMENSIONS,err,error,*999)
                        !calculate number of components with one component for each dimension and one for pressure
                        DEPENDENT_FIELD_NUMBER_OF_COMPONENTS=NUMBER_OF_DIMENSIONS+1
                        CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                          & FIELD_U_VARIABLE_TYPE,DEPENDENT_FIELD_NUMBER_OF_COMPONENTS,err,error,*999)
                        CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                          & FIELD_DELUDELN_VARIABLE_TYPE,DEPENDENT_FIELD_NUMBER_OF_COMPONENTS,err,error,*999)
                        CALL FIELD_COMPONENT_MESH_COMPONENT_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                          & 1,GEOMETRIC_MESH_COMPONENT,err,error,*999)
                        !Default to the geometric interpolation setup
                        DO I=1,DEPENDENT_FIELD_NUMBER_OF_COMPONENTS
                          CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                            & FIELD_U_VARIABLE_TYPE,I,GEOMETRIC_MESH_COMPONENT,err,error,*999)
                          CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                            & FIELD_DELUDELN_VARIABLE_TYPE,I,GEOMETRIC_MESH_COMPONENT,err,error,*999)
                        END DO
                        SELECT CASE(EQUATIONS_SET%SOLUTION_METHOD)
                          !Specify fem solution method
                          CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                            DO I=1,DEPENDENT_FIELD_NUMBER_OF_COMPONENTS
                              CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                              & FIELD_U_VARIABLE_TYPE,I,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                              CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                              & FIELD_DELUDELN_VARIABLE_TYPE,I,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                            END DO
                            CALL FIELD_SCALING_TYPE_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_SCALING_TYPE, &
                              & err,error,*999)
                            CALL FIELD_SCALING_TYPE_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,GEOMETRIC_SCALING_TYPE, &
                              & err,error,*999)
                            !Other solutions not defined yet
                          CASE DEFAULT
                            localError="The solution method of " &
                              & //TRIM(NumberToVString(EQUATIONS_SET%SOLUTION_METHOD,"*",err,error))// " is invalid."
                            CALL FlagError(localError,err,error,*999)
                        END SELECT
                      ELSE
                        !Check the user specified field
                        CALL FIELD_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_GENERAL_TYPE,err,error,*999)
                        CALL FIELD_DEPENDENT_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DEPENDENT_TYPE,err,error,*999)
                        CALL FIELD_NUMBER_OF_VARIABLES_CHECK(EQUATIONS_SET_SETUP%FIELD,2,err,error,*999)
                        CALL FIELD_VARIABLE_TYPES_CHECK(EQUATIONS_SET_SETUP%FIELD,[FIELD_U_VARIABLE_TYPE, &
                          & FIELD_DELUDELN_VARIABLE_TYPE],err,error,*999)
                        CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
                          & err,error,*999)
                        CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                          & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
                        CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
                        CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE,FIELD_DP_TYPE, &
                          & err,error,*999)
                        CALL FIELD_NUMBER_OF_COMPONENTS_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                          & NUMBER_OF_DIMENSIONS,err,error,*999)
                        !calculate number of components with one component for each dimension and one for pressure
                        DEPENDENT_FIELD_NUMBER_OF_COMPONENTS=NUMBER_OF_DIMENSIONS+1
                        CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE, &
                          & DEPENDENT_FIELD_NUMBER_OF_COMPONENTS,err,error,*999)
                        CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                          & DEPENDENT_FIELD_NUMBER_OF_COMPONENTS,err,error,*999)
                        SELECT CASE(EQUATIONS_SET%SOLUTION_METHOD)
                          CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                            CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,1, &
                              & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                            CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE,1, &
                              & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                          CASE DEFAULT
                            localError="The solution method of "//TRIM(NumberToVString(EQUATIONS_SET%SOLUTION_METHOD, &
                              &"*",err,error))//" is invalid."
                             CALL FlagError(localError,err,error,*999)
                        END SELECT
                      ENDIF
                     !Specify finish action
                    CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
                      IF(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD_AUTO_CREATED) THEN
                        CALL FIELD_CREATE_FINISH(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,err,error,*999)
                      ENDIF
                    CASE DEFAULT
                      localError="The action type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",err,error))// &
                      & " for a setup type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",err,error))// &
                      & " is invalid for a standard Stokes fluid"
                      CALL FlagError(localError,err,error,*999)
                  END SELECT
                CASE DEFAULT
                  localError="The third equations set specification of "// &
                    & TRIM(NumberToVstring(EQUATIONS_SET%SPECIFICATION(3),"*", &
                    & err,error))//" is invalid for a Stokes flow equations set."
                   CALL FLAG_ERROR(localError,err,error,*999)
                 END SELECT
            CASE(EQUATIONS_SET_SETUP_INDEPENDENT_TYPE)
              !define an independent field for ALE information
              SELECT CASE(EQUATIONS_SET%SPECIFICATION(3))
              CASE(EQUATIONS_SET_ALE_STOKES_SUBTYPE,EQUATIONS_SET_PGM_STOKES_SUBTYPE)
                SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
                 !Set start action
                  CASE(EQUATIONS_SET_SETUP_START_ACTION)
                    IF(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD_AUTO_CREATED) THEN
                      !Create the auto created independent field
                      !start field creation with name 'INDEPENDENT_FIELD'
                      CALL FIELD_CREATE_START(EQUATIONS_SET_SETUP%FIELD_USER_NUMBER,EQUATIONS_SET%REGION, &
                        & EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,err,error,*999)
                      !start creation of a new field
                      CALL FIELD_TYPE_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_GENERAL_TYPE,err,error,*999)
                      !label the field
                      CALL FIELD_LABEL_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,"Independent Field",err,error,*999)
                      !define new created field to be independent
                      CALL FIELD_DEPENDENT_TYPE_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD, &
                        & FIELD_INDEPENDENT_TYPE,err,error,*999)
                      !look for decomposition rule already defined
                      CALL FIELD_MESH_DECOMPOSITION_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_DECOMPOSITION, &
                        & err,error,*999)
                      !apply decomposition rule found on new created field
                      CALL FIELD_MESH_DECOMPOSITION_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD, &
                        & GEOMETRIC_DECOMPOSITION,err,error,*999)
                      !point new field to geometric field
                      CALL FIELD_GEOMETRIC_FIELD_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,EQUATIONS_SET% &
                        & GEOMETRY%GEOMETRIC_FIELD,err,error,*999)
                      !set number of variables to 1 (1 for U)
                      INDEPENDENT_FIELD_NUMBER_OF_VARIABLES=1
                      CALL FIELD_NUMBER_OF_VARIABLES_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD, &
                      & INDEPENDENT_FIELD_NUMBER_OF_VARIABLES,err,error,*999)
                      CALL FIELD_VARIABLE_TYPES_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD, &
                        & [FIELD_U_VARIABLE_TYPE],err,error,*999)
                      CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                        & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
                      CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                        & FIELD_DP_TYPE,err,error,*999)
                      CALL FIELD_NUMBER_OF_COMPONENTS_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                        & NUMBER_OF_DIMENSIONS,err,error,*999)
                      !calculate number of components with one component for each dimension
                      INDEPENDENT_FIELD_NUMBER_OF_COMPONENTS=NUMBER_OF_DIMENSIONS
                      CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD, &
                        & FIELD_U_VARIABLE_TYPE,INDEPENDENT_FIELD_NUMBER_OF_COMPONENTS,err,error,*999)
                      CALL FIELD_COMPONENT_MESH_COMPONENT_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                        & 1,GEOMETRIC_MESH_COMPONENT,err,error,*999)
                      !Default to the geometric interpolation setup
                      DO I=1,INDEPENDENT_FIELD_NUMBER_OF_COMPONENTS
                        CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD, &
                          & FIELD_U_VARIABLE_TYPE,I,GEOMETRIC_MESH_COMPONENT,err,error,*999)
                      END DO
                        SELECT CASE(EQUATIONS_SET%SOLUTION_METHOD)
                          !Specify fem solution method
                          CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                            DO I=1,INDEPENDENT_FIELD_NUMBER_OF_COMPONENTS
                              CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD, &
                              & FIELD_U_VARIABLE_TYPE,I,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                            END DO
                            CALL FIELD_SCALING_TYPE_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_SCALING_TYPE, &
                              & err,error,*999)
                            CALL FIELD_SCALING_TYPE_SET(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,GEOMETRIC_SCALING_TYPE, &
                              & err,error,*999)
                            !Other solutions not defined yet
                          CASE DEFAULT
                            localError="The solution method of " &
                              & //TRIM(NumberToVString(EQUATIONS_SET%SOLUTION_METHOD,"*",err,error))// " is invalid."
                            CALL FlagError(localError,err,error,*999)
                        END SELECT
                      ELSE
                        !Check the user specified field
                        CALL FIELD_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_GENERAL_TYPE,err,error,*999)
                        CALL FIELD_DEPENDENT_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_INDEPENDENT_TYPE,err,error,*999)
                        CALL FIELD_NUMBER_OF_VARIABLES_CHECK(EQUATIONS_SET_SETUP%FIELD,1,err,error,*999)
                        CALL FIELD_VARIABLE_TYPES_CHECK(EQUATIONS_SET_SETUP%FIELD,[FIELD_U_VARIABLE_TYPE],err,error,*999)
                        CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
                          & err,error,*999)
                        CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
                        CALL FIELD_NUMBER_OF_COMPONENTS_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                          & NUMBER_OF_DIMENSIONS,err,error,*999)
                        !calculate number of components with one component for each dimension and one for pressure
                        INDEPENDENT_FIELD_NUMBER_OF_COMPONENTS=NUMBER_OF_DIMENSIONS
                        CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE, &
                          & INDEPENDENT_FIELD_NUMBER_OF_COMPONENTS,err,error,*999)
                        SELECT CASE(EQUATIONS_SET%SOLUTION_METHOD)
                          CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                            CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,1, &
                              & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                            CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE,1, &
                              & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                          CASE DEFAULT
                            localError="The solution method of "//TRIM(NumberToVString(EQUATIONS_SET%SOLUTION_METHOD, &
                              &"*",err,error))//" is invalid."
                             CALL FlagError(localError,err,error,*999)
                        END SELECT
                      ENDIF
                  !Specify finish action
                  CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
                    IF(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD_AUTO_CREATED) THEN
                      CALL FIELD_CREATE_FINISH(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,err,error,*999)
                      CALL FIELD_PARAMETER_SET_CREATE(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                         & FIELD_MESH_DISPLACEMENT_SET_TYPE,err,error,*999)
                      CALL FIELD_PARAMETER_SET_CREATE(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                         & FIELD_MESH_VELOCITY_SET_TYPE,err,error,*999)
                      CALL FIELD_PARAMETER_SET_CREATE(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                        & FIELD_BOUNDARY_SET_TYPE,err,error,*999)
                    ENDIF
                    CASE DEFAULT
                      localError="The action type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",err,error))// &
                      & " for a setup type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",err,error))// &
                      & " is invalid for a standard Stokes fluid"
                      CALL FlagError(localError,err,error,*999)
                  END SELECT
                CASE DEFAULT
                  localError="The third equations set specification of "// &
                    & TRIM(NumberToVstring(EQUATIONS_SET%SPECIFICATION(3),"*", &
                    & err,error))//" is invalid for a Stokes flow equations set."
                   CALL FlagError(localError,err,error,*999)
                 END SELECT
            !Define analytic part for Stokes problem
            CASE(EQUATIONS_SET_SETUP_ANALYTIC_TYPE)
              SELECT CASE(EQUATIONS_SET%SPECIFICATION(3))
                CASE(EQUATIONS_SET_STATIC_STOKES_SUBTYPE,EQUATIONS_SET_LAPLACE_STOKES_SUBTYPE, &
                  & EQUATIONS_SET_TRANSIENT_STOKES_SUBTYPE)
                  SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
                    !Set start action
                    CASE(EQUATIONS_SET_SETUP_START_ACTION)
                      IF(EQUATIONS_SET%DEPENDENT%DEPENDENT_FINISHED) THEN
                        IF(ASSOCIATED(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD)) THEN
                          IF(ASSOCIATED(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD)) THEN
                            CALL FIELD_NUMBER_OF_COMPONENTS_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                              & NUMBER_OF_DIMENSIONS,err,error,*999)
                            SELECT CASE(EQUATIONS_SET_SETUP%ANALYTIC_FUNCTION_TYPE)
                              CASE(EQUATIONS_SET_STOKES_EQUATION_TWO_DIM_1)
                                !Set analtyic function type
                                EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE=EQUATIONS_SET_STOKES_EQUATION_TWO_DIM_1
                              CASE(EQUATIONS_SET_STOKES_EQUATION_TWO_DIM_2)
                                !Set analtyic function type
                                EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE=EQUATIONS_SET_STOKES_EQUATION_TWO_DIM_2
                              CASE(EQUATIONS_SET_STOKES_EQUATION_TWO_DIM_3)
                                !Set analtyic function type
                                EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE=EQUATIONS_SET_STOKES_EQUATION_TWO_DIM_3
                              CASE(EQUATIONS_SET_STOKES_EQUATION_TWO_DIM_4)
                                !Set analtyic function type
                                EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE=EQUATIONS_SET_STOKES_EQUATION_TWO_DIM_4
                              CASE(EQUATIONS_SET_STOKES_EQUATION_TWO_DIM_5)
                                !Set analtyic function type
                                EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE=EQUATIONS_SET_STOKES_EQUATION_TWO_DIM_5
                              CASE(EQUATIONS_SET_STOKES_EQUATION_THREE_DIM_1)
                                !Set analtyic function type
                                EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE=EQUATIONS_SET_STOKES_EQUATION_THREE_DIM_1
                              CASE(EQUATIONS_SET_STOKES_EQUATION_THREE_DIM_2)
                                !Set analtyic function type
                                EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE=EQUATIONS_SET_STOKES_EQUATION_THREE_DIM_2
                              CASE(EQUATIONS_SET_STOKES_EQUATION_THREE_DIM_3)
                                !Set analtyic function type
                                EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE=EQUATIONS_SET_STOKES_EQUATION_THREE_DIM_3
                              CASE(EQUATIONS_SET_STOKES_EQUATION_THREE_DIM_4)
                                !Set analtyic function type
                                EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE=EQUATIONS_SET_STOKES_EQUATION_THREE_DIM_4
                              CASE(EQUATIONS_SET_STOKES_EQUATION_THREE_DIM_5)
                                !Set analtyic function type
                                EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE=EQUATIONS_SET_STOKES_EQUATION_THREE_DIM_5
                              CASE DEFAULT
                                localError="The specified analytic function type of "// &
                                  & TRIM(NumberToVString(EQUATIONS_SET_SETUP%ANALYTIC_FUNCTION_TYPE,"*",err,error))// &
                                  & " is invalid for an analytic Stokes problem."
                                CALL FlagError(localError,err,error,*999)
                            END SELECT
                          ELSE
                            CALL FlagError("Equations set geometric field is not associated.",err,error,*999)
                          ENDIF
                        ELSE
                          CALL FlagError("Equations set dependent field is not associated.",err,error,*999)
                        ENDIF
                      ELSE
                        CALL FlagError("Equations set dependent field has not been finished.",err,error,*999)
                      ENDIF
                    CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
                      IF(ASSOCIATED(EQUATIONS_SET%ANALYTIC)) THEN
                        IF(ASSOCIATED(EQUATIONS_SET%ANALYTIC%ANALYTIC_FIELD)) THEN
                          IF(EQUATIONS_SET%ANALYTIC%ANALYTIC_FIELD_AUTO_CREATED) THEN
                            CALL FIELD_CREATE_FINISH(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,err,error,*999)
                          ENDIF
                        ENDIF
                      ELSE
                        CALL FlagError("Equations set analytic is not associated.",err,error,*999)
                      ENDIF
                    CASE DEFAULT
                      localError="The action type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",err,error))// &
                        & " for a setup type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",err,error))// &
                        & " is invalid for an analytic Stokes problem."
                      CALL FlagError(localError,err,error,*999)
                    END SELECT
                  CASE DEFAULT
                    localError="The third equations set specification of "// &
                      & TRIM(NumberToVString(EQUATIONS_SET%SPECIFICATION(3),"*",err,error))// &
                       " is invalid for a Stokes flow equations set."
                     CALL FLAG_ERROR(localError,err,error,*999)
                END SELECT
            !Define materials field
            CASE(EQUATIONS_SET_SETUP_MATERIALS_TYPE)
              SELECT CASE(EQUATIONS_SET%SPECIFICATION(3))
                CASE(EQUATIONS_SET_STATIC_STOKES_SUBTYPE,EQUATIONS_SET_LAPLACE_STOKES_SUBTYPE, &
                & EQUATIONS_SET_TRANSIENT_STOKES_SUBTYPE,EQUATIONS_SET_ALE_STOKES_SUBTYPE, &
                & EQUATIONS_SET_PGM_STOKES_SUBTYPE)
                  !variable X with has Y components, here Y represents viscosity only
                  MATERIAL_FIELD_NUMBER_OF_VARIABLES=1!X
                  MATERIAL_FIELD_NUMBER_OF_COMPONENTS=2!Y
                  SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
                  !Specify start action
                    CASE(EQUATIONS_SET_SETUP_START_ACTION)
                      EQUATIONS_MATERIALS=>EQUATIONS_SET%MATERIALS
                      IF(ASSOCIATED(EQUATIONS_MATERIALS)) THEN
                        IF(EQUATIONS_MATERIALS%MATERIALS_FIELD_AUTO_CREATED) THEN
                          !Create the auto created materials field
                          !start field creation with name 'MATERIAL_FIELD'
                          CALL FIELD_CREATE_START(EQUATIONS_SET_SETUP%FIELD_USER_NUMBER,EQUATIONS_SET%REGION,EQUATIONS_SET% &
                            & MATERIALS%MATERIALS_FIELD,err,error,*999)
                          CALL FIELD_TYPE_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_MATERIAL_TYPE,err,error,*999)
                          !label the field
                          CALL FIELD_LABEL_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,"Materials Field",err,error,*999)
                          CALL FIELD_DEPENDENT_TYPE_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_INDEPENDENT_TYPE, &
                            & err,error,*999)
                          CALL FIELD_MESH_DECOMPOSITION_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_DECOMPOSITION, &
                            & err,error,*999)
                          !apply decomposition rule found on new created field
                          CALL FIELD_MESH_DECOMPOSITION_SET_AND_LOCK(EQUATIONS_SET%MATERIALS%MATERIALS_FIELD, &
                            & GEOMETRIC_DECOMPOSITION,err,error,*999)
                          !point new field to geometric field
                          CALL FIELD_GEOMETRIC_FIELD_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,EQUATIONS_SET%GEOMETRY% &
                            & GEOMETRIC_FIELD,err,error,*999)
                          CALL FIELD_NUMBER_OF_VARIABLES_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD, &
                            & MATERIAL_FIELD_NUMBER_OF_VARIABLES,err,error,*999)
                          CALL FIELD_VARIABLE_TYPES_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,[FIELD_U_VARIABLE_TYPE], &
                            & err,error,*999)
                          CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                            & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
                          CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                            & FIELD_DP_TYPE,err,error,*999)
                          CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                            & MATERIAL_FIELD_NUMBER_OF_COMPONENTS,err,error,*999)
                          CALL FIELD_COMPONENT_MESH_COMPONENT_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                            & 1,GEOMETRIC_COMPONENT_NUMBER,err,error,*999)
                          CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                            & 1,GEOMETRIC_COMPONENT_NUMBER,err,error,*999)
                          CALL FIELD_COMPONENT_INTERPOLATION_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                            & 1,FIELD_CONSTANT_INTERPOLATION,err,error,*999)
                        !Default the field scaling to that of the geometric field
                          CALL FIELD_SCALING_TYPE_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_SCALING_TYPE,err,error,*999)
                          CALL FIELD_SCALING_TYPE_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,GEOMETRIC_SCALING_TYPE,err,error,*999)
                        ELSE
                          !Check the user specified field
                          CALL FIELD_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_MATERIAL_TYPE,err,error,*999)
                          CALL FIELD_DEPENDENT_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_INDEPENDENT_TYPE,err,error,*999)
                          CALL FIELD_NUMBER_OF_VARIABLES_CHECK(EQUATIONS_SET_SETUP%FIELD,1,err,error,*999)
                          CALL FIELD_VARIABLE_TYPES_CHECK(EQUATIONS_SET_SETUP%FIELD,[FIELD_U_VARIABLE_TYPE],err,error,*999)
                          CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
                            & err,error,*999)
                          CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
                          CALL FIELD_NUMBER_OF_COMPONENTS_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                            & NUMBER_OF_DIMENSIONS,err,error,*999)
                          CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,1,err,error,*999)
                        ENDIF
                      ELSE
                        CALL FlagError("Equations set materials is not associated.",err,error,*999)
                      END IF
                    !Specify start action
                    CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
                      EQUATIONS_MATERIALS=>EQUATIONS_SET%MATERIALS
                      IF(ASSOCIATED(EQUATIONS_MATERIALS)) THEN
                        IF(EQUATIONS_MATERIALS%MATERIALS_FIELD_AUTO_CREATED) THEN
                          !Finish creating the materials field
                          CALL FIELD_CREATE_FINISH(EQUATIONS_MATERIALS%MATERIALS_FIELD,err,error,*999)
                          !Set the default values for the materials field
                          !First set the mu values to 0.001
                          !MATERIAL_FIELD_NUMBER_OF_COMPONENTS
                          ! viscosity=1
                          CALL FIELD_COMPONENT_VALUES_INITIALISE(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                            & FIELD_VALUES_SET_TYPE,1,1.0_DP,err,error,*999)
                          ! density=2
!\todo: Initialise fields properly
  !                       CALL FIELD_COMPONENT_VALUES_INITIALISE(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
  !                         & FIELD_VALUES_SET_TYPE,2,100.0_DP,err,error,*999)
                        ENDIF
                      ELSE
                        CALL FlagError("Equations set materials is not associated.",err,error,*999)
                      ENDIF
                    CASE DEFAULT
                      localError="The action type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",err,error))// &
                      & " for a setup type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",err,error))// &
                        & " is invalid for Stokes equation."
                      CALL FlagError(localError,err,error,*999)
                  END SELECT
                CASE DEFAULT
                  localError="The third equations set specification of "// &
                    & TRIM(NumberToVString(EQUATIONS_SET%SPECIFICATION(3),"*",err,error))// &
                    & " is invalid for a Stokes flow equations set."
                  CALL FlagError(localError,err,error,*999)
              END SELECT
            !Define the source field
            CASE(EQUATIONS_SET_SETUP_SOURCE_TYPE)
              SELECT CASE(EQUATIONS_SET%SPECIFICATION(3))
                CASE(EQUATIONS_SET_STATIC_STOKES_SUBTYPE,EQUATIONS_SET_LAPLACE_STOKES_SUBTYPE, &
                  & EQUATIONS_SET_TRANSIENT_STOKES_SUBTYPE,EQUATIONS_SET_ALE_STOKES_SUBTYPE, &
                  & EQUATIONS_SET_PGM_STOKES_SUBTYPE)
                  !TO DO: INCLUDE GRAVITY AS SOURCE TYPE
                  SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
                  CASE(EQUATIONS_SET_SETUP_START_ACTION)
                      !Do nothing
                    CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
                      !Do nothing
                      !? Maybe set finished flag????
                    CASE DEFAULT
                      localError="The action type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",err,error))// &
                        & " for a setup type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",err,error))// &
                        & " is invalid for a standard Stokes fluid."
                      CALL FlagError(localError,err,error,*999)
                  END SELECT
                CASE DEFAULT
                  localError="The third equations set specification of "// &
                    & TRIM(NumberToVString(EQUATIONS_SET%SPECIFICATION(3),"*",err,error))// &
                    & " is invalid for a Stokes flow equations set."
                  CALL FLAG_ERROR(localError,err,error,*999)
              END SELECT
            !Define equations type
            CASE(EQUATIONS_SET_SETUP_EQUATIONS_TYPE)
              SELECT CASE(EQUATIONS_SET%SPECIFICATION(3))
                CASE(EQUATIONS_SET_STATIC_STOKES_SUBTYPE,EQUATIONS_SET_LAPLACE_STOKES_SUBTYPE)
                  SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
                    CASE(EQUATIONS_SET_SETUP_START_ACTION)
                      EQUATIONS_MATERIALS=>EQUATIONS_SET%MATERIALS
                      IF(ASSOCIATED(EQUATIONS_MATERIALS)) THEN
                        IF(EQUATIONS_MATERIALS%MATERIALS_FINISHED) THEN
                          CALL Equations_CreateStart(EQUATIONS_SET,equations,err,error,*999)
                          CALL Equations_LinearityTypeSet(equations,EQUATIONS_LINEAR,err,error,*999)
                          CALL Equations_TimeDependenceTypeSet(equations,EQUATIONS_STATIC,err,error,*999)
                        ELSE
                          CALL FlagError("Equations set materials has not been finished.",err,error,*999)
                        ENDIF
                      ELSE
                        CALL FlagError("Equations materials is not associated.",err,error,*999)
                      ENDIF
                    CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
                      SELECT CASE(EQUATIONS_SET%SOLUTION_METHOD)
                        CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                        !Finish the equations creation
                          CALL EquationsSet_EquationsGet(EQUATIONS_SET,equations,err,error,*999)
                          CALL Equations_CreateFinish(equations,err,error,*999)
                          NULLIFY(vectorEquations)
                          CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
                          !Create the equations mapping.
                          CALL EquationsMapping_VectorCreateStart(vectorEquations,FIELD_DELUDELN_VARIABLE_TYPE,vectorMapping, &
                            & err,error,*999)
                          CALL EquationsMapping_LinearMatricesNumberSet(vectorMapping,1,err,error,*999)
                          CALL EquationsMapping_LinearMatricesVariableTypesSet(vectorMapping,[FIELD_U_VARIABLE_TYPE], &
                            & err,error,*999)
                          CALL EquationsMapping_RHSVariableTypeSet(vectorMapping,FIELD_DELUDELN_VARIABLE_TYPE, &
                            & err,error,*999)
                          CALL EquationsMapping_VectorCreateFinish(vectorMapping,err,error,*999)
                          !Create the equations matrices
                          CALL EquationsMatrices_VectorCreateStart(vectorEquations,vectorMatrices,err,error,*999)
                          SELECT CASE(equations%sparsityType)
                            CASE(EQUATIONS_MATRICES_FULL_MATRICES)
                              CALL EquationsMatrices_LinearStorageTypeSet(vectorMatrices,[MATRIX_BLOCK_STORAGE_TYPE], &
                                & err,error,*999)
                            CASE(EQUATIONS_MATRICES_SPARSE_MATRICES)
                              CALL EquationsMatrices_LinearStorageTypeSet(vectorMatrices, &
                                & [MATRIX_COMPRESSED_ROW_STORAGE_TYPE],err,error,*999)
                              CALL EquationsMatrices_LinearStructureTypeSet(vectorMatrices, &
                              & [EQUATIONS_MATRIX_FEM_STRUCTURE],err,error,*999)
                            CASE DEFAULT
                              localError="The equations matrices sparsity type of "// &
                                & TRIM(NumberToVString(equations%sparsityType,"*",err,error))//" is invalid."
                              CALL FlagError(localError,err,error,*999)
                          END SELECT
                          CALL EquationsMatrices_VectorCreateFinish(vectorMatrices,err,error,*999)
                        CASE DEFAULT
                          localError="The solution method of "//TRIM(NumberToVString(EQUATIONS_SET%SOLUTION_METHOD, &
                            & "*",err,error))//" is invalid."
                          CALL FlagError(localError,err,error,*999)
                      END SELECT
                    CASE DEFAULT
                      localError="The action type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",err,error))// &
                        & " for a setup type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",err,error))// &
                        & " is invalid for a steady Laplace equation."
                      CALL FlagError(localError,err,error,*999)
                  END SELECT
                CASE(EQUATIONS_SET_TRANSIENT_STOKES_SUBTYPE,EQUATIONS_SET_ALE_STOKES_SUBTYPE,EQUATIONS_SET_PGM_STOKES_SUBTYPE)
                SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
                    CASE(EQUATIONS_SET_SETUP_START_ACTION)
                      EQUATIONS_MATERIALS=>EQUATIONS_SET%MATERIALS
                      IF(ASSOCIATED(EQUATIONS_MATERIALS)) THEN
                        IF(EQUATIONS_MATERIALS%MATERIALS_FINISHED) THEN
                          CALL Equations_CreateStart(EQUATIONS_SET,equations,err,error,*999)
                          CALL Equations_EquationTypeSet(equations,EQUATIONS_VECTOR_TYPE,err,error,*999)
                          CALL Equations_LinearityTypeSet(equations,EQUATIONS_LINEAR,err,error,*999)
                          CALL Equations_TimeDependenceTypeSet(equations,EQUATIONS_FIRST_ORDER_DYNAMIC,err,error,*999)
                        ELSE
                          CALL FlagError("Equations set materials has not been finished.",err,error,*999)
                        ENDIF
                      ELSE
                        CALL FlagError("Equations materials is not associated.",err,error,*999)
                      ENDIF
                    CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
                      SELECT CASE(EQUATIONS_SET%SOLUTION_METHOD)
                        CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                          !Finish the equations creation
                          CALL EquationsSet_EquationsGet(EQUATIONS_SET,equations,err,error,*999)
                          CALL Equations_CreateFinish(equations,err,error,*999)
                          NULLIFY(vectorEquations)
                          CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
                          !Create the equations mapping.
                          CALL EquationsMapping_VectorCreateStart(vectorEquations,FIELD_DELUDELN_VARIABLE_TYPE,vectorMapping, &
                            & err,error,*999)
                          CALL EquationsMapping_DynamicMatricesSet(vectorMapping,.TRUE.,.TRUE.,err,error,*999)
                          CALL EquationsMapping_DynamicVariableTypeSet(vectorMapping,FIELD_U_VARIABLE_TYPE,err,error,*999)
                          CALL EquationsMapping_RHSVariableTypeSet(vectorMapping,FIELD_DELUDELN_VARIABLE_TYPE, &
                            & err,error,*999)
                          CALL EquationsMapping_VectorCreateFinish(vectorMapping,err,error,*999)
                          !Create the equations matrices
                          CALL EquationsMatrices_VectorCreateStart(vectorEquations,vectorMatrices,err,error,*999)
                          !Set up matrix storage and structure
                          IF(equations%lumpingType==EQUATIONS_LUMPED_MATRICES) THEN
                            !Set up lumping
                            CALL EquationsMatrices_DynamicLumpingTypeSet(vectorMatrices, &
                              & [EQUATIONS_MATRIX_UNLUMPED,EQUATIONS_MATRIX_LUMPED],err,error,*999)
                            CALL EquationsMatrices_DynamicStorageTypeSet(vectorMatrices, &
                              & [DISTRIBUTED_MATRIX_COMPRESSED_ROW_STORAGE_TYPE,DISTRIBUTED_MATRIX_DIAGONAL_STORAGE_TYPE] &
                              & ,err,error,*999)
                            CALL EquationsMatrices_DynamicStructureTypeSet(vectorMatrices, &
                              & [EQUATIONS_MATRIX_FEM_STRUCTURE,EQUATIONS_MATRIX_DIAGONAL_STRUCTURE],err,error,*999)
                          ELSE
                            SELECT CASE(equations%sparsityType)
                              CASE(EQUATIONS_MATRICES_FULL_MATRICES)
                              CALL EquationsMatrices_LinearStorageTypeSet(vectorMatrices, &
                                  & [DISTRIBUTED_MATRIX_BLOCK_STORAGE_TYPE,DISTRIBUTED_MATRIX_BLOCK_STORAGE_TYPE],err,error,*999)
                              CASE(EQUATIONS_MATRICES_SPARSE_MATRICES)
                                CALL EquationsMatrices_DynamicStorageTypeSet(vectorMatrices, &
                                  & [DISTRIBUTED_MATRIX_COMPRESSED_ROW_STORAGE_TYPE, &
                                  & DISTRIBUTED_MATRIX_COMPRESSED_ROW_STORAGE_TYPE],err,error,*999)
                                CALL EquationsMatrices_DynamicStructureTypeSet(vectorMatrices, &
                                  & [EQUATIONS_MATRIX_FEM_STRUCTURE,EQUATIONS_MATRIX_FEM_STRUCTURE],err,error,*999)
                              CASE DEFAULT
                                localError="The equations matrices sparsity type of "// &
                                  & TRIM(NumberToVString(equations%sparsityType,"*",err,error))//" is invalid."
                                CALL FlagError(localError,err,error,*999)
                            END SELECT
                          ENDIF
                          CALL EquationsMatrices_VectorCreateFinish(vectorMatrices,err,error,*999)
                        CASE DEFAULT
                          localError="The solution method of "//TRIM(NumberToVString(EQUATIONS_SET%SOLUTION_METHOD,"*", &
                            & err,error))//" is invalid."
                          CALL FlagError(localError,err,error,*999)
                      END SELECT
                    CASE DEFAULT
                    localError="The action type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",err,error))// &
                        & " for a setup type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",err,error))// &
                        & " is invalid for a Stokes equation."
                      CALL FlagError(localError,err,error,*999)
                  END SELECT
                CASE DEFAULT
                  localError="The third equations set specification of "// &
                    & TRIM(NumberToVString(EQUATIONS_SET%SPECIFICATION(3),"*",err,error))// &
                    & " is invalid for a Stokes flow equations set."
                  CALL FlagError(localError,err,error,*999)
              END SELECT
            CASE DEFAULT
              localError="The setup type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",err,error))// &
                & " is invalid for a standard Stokes fluid."
              CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE DEFAULT
          localError="The equations set subtype of "//TRIM(NumberToVString(EQUATIONS_SET%SPECIFICATION(3),"*",err,error))// &
            & " does not equal a standard Stokes fluid subtype."
          CALL FlagError(localError,err,error,*999)
      END SELECT
    ELSE
      CALL FlagError("Equations set is not associated.",err,error,*999)
    ENDIF

    EXITS("STOKES_EQUATIONS_SET_SETUP")
    RETURN
999 ERRORSEXITS("STOKES_EQUATIONS_SET_SETUP",err,error)
    RETURN 1
  END SUBROUTINE STOKES_EQUATIONS_SET_SETUP

!
!================================================================================================================================
!

  !>Sets the problem specification for a Stokes fluid problem.
  SUBROUTINE Stokes_ProblemSpecificationSet(problem,problemSpecification,err,error,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: problem !<A pointer to the problem to set the problem specification for
    INTEGER(INTG), INTENT(IN) :: problemSpecification(:) !<The problem specification to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
    INTEGER(INTG) :: problemSubtype

    ENTERS("Stokes_ProblemSpecificationSet",err,error,*999)

    IF(ASSOCIATED(problem)) THEN
      IF(SIZE(problemSpecification,1)==3) THEN
        problemSubtype=problemSpecification(3)
        SELECT CASE(problemSubtype)
        CASE(PROBLEM_STATIC_STOKES_SUBTYPE, &
            & PROBLEM_LAPLACE_STOKES_SUBTYPE, &
            & PROBLEM_TRANSIENT_STOKES_SUBTYPE, &
            & PROBLEM_ALE_STOKES_SUBTYPE, &
            & PROBLEM_PGM_STOKES_SUBTYPE)
          !All ok
        CASE(PROBLEM_OPTIMISED_STOKES_SUBTYPE)
          CALL FlagError("Not implemented yet.",err,error,*999)
        CASE DEFAULT
          localError="The third problem specification of "//TRIM(NumberToVstring(problemSubtype,"*",err,error))// &
            & " is not valid for a Stokes flow fluid mechanics problem."
          CALL FlagError(localError,err,error,*999)
        END SELECT
        IF(ALLOCATED(problem%specification)) THEN
          CALL FlagError("Problem specification is already allocated.",err,error,*999)
        ELSE
          ALLOCATE(problem%specification(3),stat=err)
          IF(err/=0) CALL FlagError("Could not allocate problem specification.",err,error,*999)
        ENDIF
        problem%specification(1:3)=[PROBLEM_FLUID_MECHANICS_CLASS,PROBLEM_STOKES_EQUATION_TYPE,problemSubtype]
      ELSE
        CALL FlagError("Stokes flow problem specification must have >= 3 entries.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Problem is not associated.",err,error,*999)
    ENDIF

    EXITS("Stokes_ProblemSpecificationSet")
    RETURN
999 ERRORSEXITS("Stokes_ProblemSpecificationSet",err,error)
    RETURN 1

  END SUBROUTINE Stokes_ProblemSpecificationSet

!
!================================================================================================================================
!

  !>Sets up the Stokes problem.
  SUBROUTINE STOKES_PROBLEM_SETUP(PROBLEM,PROBLEM_SETUP,err,error,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM !<A pointer to the problem set to setup a Stokes fluid on.
    TYPE(PROBLEM_SETUP_TYPE), INTENT(INOUT) :: PROBLEM_SETUP !<The problem setup information
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP,CONTROL_LOOP_ROOT
    TYPE(SOLVER_TYPE), POINTER :: SOLVER, MESH_SOLVER
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS,MESH_SOLVER_EQUATIONS
    TYPE(SOLVERS_TYPE), POINTER :: SOLVERS

    ENTERS("STOKES_PROBLEM_SETUP",err,error,*999)

    NULLIFY(CONTROL_LOOP)
    NULLIFY(SOLVER)
    NULLIFY(MESH_SOLVER)
    NULLIFY(SOLVER_EQUATIONS)
    NULLIFY(MESH_SOLVER_EQUATIONS)
    NULLIFY(SOLVERS)
    IF(ASSOCIATED(PROBLEM)) THEN
      IF(.NOT.ALLOCATED(problem%specification)) THEN
        CALL FlagError("Problem specification is not allocated.",err,error,*999)
      ELSE IF(SIZE(problem%specification,1)<3) THEN
        CALL FlagError("Problem specification must have three entries for a Stokes problem.",err,error,*999)
      END IF
      SELECT CASE(PROBLEM%SPECIFICATION(3))
        !Set problem subtype for steady state Stokes problems
        CASE(PROBLEM_STATIC_STOKES_SUBTYPE,PROBLEM_LAPLACE_STOKES_SUBTYPE)
          SELECT CASE(PROBLEM_SETUP%SETUP_TYPE)
            CASE(PROBLEM_SETUP_INITIAL_TYPE)
              SELECT CASE(PROBLEM_SETUP%ACTION_TYPE)
                CASE(PROBLEM_SETUP_START_ACTION)
                  !Do nothing????
                CASE(PROBLEM_SETUP_FINISH_ACTION)
                  !Do nothing???
                CASE DEFAULT
                  localError="The action type of "//TRIM(NumberToVString(PROBLEM_SETUP%ACTION_TYPE,"*",err,error))// &
                    & " for a setup type of "//TRIM(NumberToVString(PROBLEM_SETUP%SETUP_TYPE,"*",err,error))// &
                    & " is invalid for a standard Stokes fluid."
                  CALL FlagError(localError,err,error,*999)
              END SELECT
            CASE(PROBLEM_SETUP_CONTROL_TYPE)
              SELECT CASE(PROBLEM_SETUP%ACTION_TYPE)
                CASE(PROBLEM_SETUP_START_ACTION)
                  !Set up a simple control loop
                  CALL CONTROL_LOOP_CREATE_START(PROBLEM,CONTROL_LOOP,err,error,*999)
                CASE(PROBLEM_SETUP_FINISH_ACTION)
                  !Finish the control loops
                  CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
                  CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,err,error,*999)
                  CALL CONTROL_LOOP_CREATE_FINISH(CONTROL_LOOP,err,error,*999)
                CASE DEFAULT
                  localError="The action type of "//TRIM(NumberToVString(PROBLEM_SETUP%ACTION_TYPE,"*",err,error))// &
                    & " for a setup type of "//TRIM(NumberToVString(PROBLEM_SETUP%SETUP_TYPE,"*",err,error))// &
                    & " is invalid for a standard Stokes fluid."
                  CALL FlagError(localError,err,error,*999)
              END SELECT
            CASE(PROBLEM_SETUP_SOLVERS_TYPE)
              !Get the control loop
              CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
              CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,err,error,*999)
              SELECT CASE(PROBLEM_SETUP%ACTION_TYPE)
                CASE(PROBLEM_SETUP_START_ACTION)
                  !Start the solvers creation
                  CALL SOLVERS_CREATE_START(CONTROL_LOOP,SOLVERS,err,error,*999)
                  CALL SOLVERS_NUMBER_SET(SOLVERS,1,err,error,*999)
                  !Set the solver to be a linear solver
                  CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,err,error,*999)
                  CALL SOLVER_TYPE_SET(SOLVER,SOLVER_LINEAR_TYPE,err,error,*999)
                  !Set solver defaults
                  CALL SOLVER_LIBRARY_TYPE_SET(SOLVER,SOLVER_PETSC_LIBRARY,err,error,*999)
                CASE(PROBLEM_SETUP_FINISH_ACTION)
                  !Get the solvers
                  CALL CONTROL_LOOP_SOLVERS_GET(CONTROL_LOOP,SOLVERS,err,error,*999)
                  !Finish the solvers creation
                  CALL SOLVERS_CREATE_FINISH(SOLVERS,err,error,*999)
                CASE DEFAULT
                  localError="The action type of "//TRIM(NumberToVString(PROBLEM_SETUP%ACTION_TYPE,"*",err,error))// &
                    & " for a setup type of "//TRIM(NumberToVString(PROBLEM_SETUP%SETUP_TYPE,"*",err,error))// &
                    & " is invalid for a standard Stokes fluid."
                  CALL FlagError(localError,err,error,*999)
              END SELECT
            CASE(PROBLEM_SETUP_SOLVER_EQUATIONS_TYPE)
              SELECT CASE(PROBLEM_SETUP%ACTION_TYPE)
                CASE(PROBLEM_SETUP_START_ACTION)
                  !Get the control loop
                  CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
                  CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,err,error,*999)
                  !Get the solver
                  CALL CONTROL_LOOP_SOLVERS_GET(CONTROL_LOOP,SOLVERS,err,error,*999)
                  CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,err,error,*999)
                  !Create the solver equations
                  CALL SOLVER_EQUATIONS_CREATE_START(SOLVER,SOLVER_EQUATIONS,err,error,*999)
                  CALL SOLVER_EQUATIONS_LINEARITY_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_LINEAR,err,error,*999)
                  CALL SOLVER_EQUATIONS_TIME_DEPENDENCE_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_STATIC,err,error,*999)
                  CALL SOLVER_EQUATIONS_SPARSITY_TYPE_SET(SOLVER_EQUATIONS,SOLVER_SPARSE_MATRICES,err,error,*999)
                CASE(PROBLEM_SETUP_FINISH_ACTION)
                  !Get the control loop
                  CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
                  CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,err,error,*999)
                  !Get the solver equations
                  CALL CONTROL_LOOP_SOLVERS_GET(CONTROL_LOOP,SOLVERS,err,error,*999)
                  CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,err,error,*999)
                  CALL SOLVER_SOLVER_EQUATIONS_GET(SOLVER,SOLVER_EQUATIONS,err,error,*999)
                  !Finish the solver equations creation
                  CALL SOLVER_EQUATIONS_CREATE_FINISH(SOLVER_EQUATIONS,err,error,*999)
                CASE DEFAULT
                  localError="The action type of "//TRIM(NumberToVString(PROBLEM_SETUP%ACTION_TYPE,"*",err,error))// &
                    & " for a setup type of "//TRIM(NumberToVString(PROBLEM_SETUP%SETUP_TYPE,"*",err,error))// &
                    & " is invalid for a standard Stokes fluid."
                  CALL FlagError(localError,err,error,*999)
              END SELECT
            CASE DEFAULT
              localError="The setup type of "//TRIM(NumberToVString(PROBLEM_SETUP%SETUP_TYPE,"*",err,error))// &
                & " is invalid for a standard Stokes fluid."
              CALL FlagError(localError,err,error,*999)
          END SELECT
        !Set problem subtype for transient Stokes problems
        CASE(PROBLEM_TRANSIENT_STOKES_SUBTYPE,PROBLEM_PGM_STOKES_SUBTYPE)
          SELECT CASE(PROBLEM_SETUP%SETUP_TYPE)
            CASE(PROBLEM_SETUP_INITIAL_TYPE)
              SELECT CASE(PROBLEM_SETUP%ACTION_TYPE)
                CASE(PROBLEM_SETUP_START_ACTION)
                  !Do nothing????
                CASE(PROBLEM_SETUP_FINISH_ACTION)
                  !Do nothing????
                CASE DEFAULT
                  localError="The action type of "//TRIM(NumberToVString(PROBLEM_SETUP%ACTION_TYPE,"*",err,error))// &
                    & " for a setup type of "//TRIM(NumberToVString(PROBLEM_SETUP%SETUP_TYPE,"*",err,error))// &
                    & " is invalid for a transient Stokes fluid."
                  CALL FlagError(localError,err,error,*999)
              END SELECT
            CASE(PROBLEM_SETUP_CONTROL_TYPE)
              SELECT CASE(PROBLEM_SETUP%ACTION_TYPE)
                CASE(PROBLEM_SETUP_START_ACTION)
                  !Set up a time control loop
                  CALL CONTROL_LOOP_CREATE_START(PROBLEM,CONTROL_LOOP,err,error,*999)
                  CALL CONTROL_LOOP_TYPE_SET(CONTROL_LOOP,PROBLEM_CONTROL_TIME_LOOP_TYPE,err,error,*999)
                CASE(PROBLEM_SETUP_FINISH_ACTION)
                  !Finish the control loops
                  CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
                  CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,err,error,*999)
                  CALL CONTROL_LOOP_CREATE_FINISH(CONTROL_LOOP,err,error,*999)
                CASE DEFAULT
                  localError="The action type of "//TRIM(NumberToVString(PROBLEM_SETUP%ACTION_TYPE,"*",err,error))// &
                    & " for a setup type of "//TRIM(NumberToVString(PROBLEM_SETUP%SETUP_TYPE,"*",err,error))// &
                    & " is invalid for a transient Stokes fluid."
                  CALL FlagError(localError,err,error,*999)
              END SELECT
            CASE(PROBLEM_SETUP_SOLVERS_TYPE)
              !Get the control loop
              CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
              CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,err,error,*999)
              SELECT CASE(PROBLEM_SETUP%ACTION_TYPE)
                CASE(PROBLEM_SETUP_START_ACTION)
                  !Start the solvers creation
                  CALL SOLVERS_CREATE_START(CONTROL_LOOP,SOLVERS,err,error,*999)
                  CALL SOLVERS_NUMBER_SET(SOLVERS,1,err,error,*999)
                  !Set the solver to be a first order dynamic solver
                  CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,err,error,*999)
                  CALL SOLVER_TYPE_SET(SOLVER,SOLVER_DYNAMIC_TYPE,err,error,*999)
                  CALL SOLVER_DYNAMIC_ORDER_SET(SOLVER,SOLVER_DYNAMIC_FIRST_ORDER,err,error,*999)
                  !Set solver defaults
                  CALL SOLVER_DYNAMIC_DEGREE_SET(SOLVER,SOLVER_DYNAMIC_FIRST_DEGREE,err,error,*999)
                  CALL SOLVER_DYNAMIC_SCHEME_SET(SOLVER,SOLVER_DYNAMIC_CRANK_NICOLSON_SCHEME,err,error,*999)
                  CALL SOLVER_LIBRARY_TYPE_SET(SOLVER,SOLVER_CMISS_LIBRARY,err,error,*999)
                CASE(PROBLEM_SETUP_FINISH_ACTION)
                  !Get the solvers
                  CALL CONTROL_LOOP_SOLVERS_GET(CONTROL_LOOP,SOLVERS,err,error,*999)
                  !Finish the solvers creation
                  CALL SOLVERS_CREATE_FINISH(SOLVERS,err,error,*999)
                CASE DEFAULT
                  localError="The action type of "//TRIM(NumberToVString(PROBLEM_SETUP%ACTION_TYPE,"*",err,error))// &
                   & " for a setup type of "//TRIM(NumberToVString(PROBLEM_SETUP%SETUP_TYPE,"*",err,error))// &
                   & " is invalid for a transient Stokes fluid."
                  CALL FlagError(localError,err,error,*999)
              END SELECT
            CASE(PROBLEM_SETUP_SOLVER_EQUATIONS_TYPE)
              SELECT CASE(PROBLEM_SETUP%ACTION_TYPE)
                CASE(PROBLEM_SETUP_START_ACTION)
                  !Get the control loop
                  CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
                  CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,err,error,*999)
                  !Get the solver
                  CALL CONTROL_LOOP_SOLVERS_GET(CONTROL_LOOP,SOLVERS,err,error,*999)
                  CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,err,error,*999)
                  !Create the solver equations
                  CALL SOLVER_EQUATIONS_CREATE_START(SOLVER,SOLVER_EQUATIONS,err,error,*999)
                  CALL SOLVER_EQUATIONS_LINEARITY_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_LINEAR,err,error,*999)
                  CALL SOLVER_EQUATIONS_TIME_DEPENDENCE_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_FIRST_ORDER_DYNAMIC,&
                  & err,error,*999)
                  CALL SOLVER_EQUATIONS_SPARSITY_TYPE_SET(SOLVER_EQUATIONS,SOLVER_SPARSE_MATRICES,err,error,*999)
                CASE(PROBLEM_SETUP_FINISH_ACTION)
                  !Get the control loop
                  CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
                  CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,err,error,*999)
                  !Get the solver equations
                  CALL CONTROL_LOOP_SOLVERS_GET(CONTROL_LOOP,SOLVERS,err,error,*999)
                  CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,err,error,*999)
                  CALL SOLVER_SOLVER_EQUATIONS_GET(SOLVER,SOLVER_EQUATIONS,err,error,*999)
                  !Finish the solver equations creation
                  CALL SOLVER_EQUATIONS_CREATE_FINISH(SOLVER_EQUATIONS,err,error,*999)
                CASE DEFAULT
                  localError="The action type of "//TRIM(NumberToVString(PROBLEM_SETUP%ACTION_TYPE,"*",err,error))// &
                    & " for a setup type of "//TRIM(NumberToVString(PROBLEM_SETUP%SETUP_TYPE,"*",err,error))// &
                    & " is invalid for a transient Stokes fluid."
                  CALL FlagError(localError,err,error,*999)
              END SELECT
            CASE DEFAULT
              localError="The setup type of "//TRIM(NumberToVString(PROBLEM_SETUP%SETUP_TYPE,"*",err,error))// &
                & " is invalid for a transient Stokes fluid."
              CALL FlagError(localError,err,error,*999)
          END SELECT
        !Set problem subtype for ALE Stokes problems
        CASE(PROBLEM_ALE_STOKES_SUBTYPE)
          SELECT CASE(PROBLEM_SETUP%SETUP_TYPE)
            CASE(PROBLEM_SETUP_INITIAL_TYPE)
              SELECT CASE(PROBLEM_SETUP%ACTION_TYPE)
                CASE(PROBLEM_SETUP_START_ACTION)
                  !Do nothing????
                CASE(PROBLEM_SETUP_FINISH_ACTION)
                  !Do nothing????
                CASE DEFAULT
                  localError="The action type of "//TRIM(NumberToVString(PROBLEM_SETUP%ACTION_TYPE,"*",err,error))// &
                    & " for a setup type of "//TRIM(NumberToVString(PROBLEM_SETUP%SETUP_TYPE,"*",err,error))// &
                    & " is invalid for a ALE Stokes fluid."
                  CALL FlagError(localError,err,error,*999)
              END SELECT
            CASE(PROBLEM_SETUP_CONTROL_TYPE)
              SELECT CASE(PROBLEM_SETUP%ACTION_TYPE)
                CASE(PROBLEM_SETUP_START_ACTION)
                  !Set up a time control loop
                  CALL CONTROL_LOOP_CREATE_START(PROBLEM,CONTROL_LOOP,err,error,*999)
                  CALL CONTROL_LOOP_TYPE_SET(CONTROL_LOOP,PROBLEM_CONTROL_TIME_LOOP_TYPE,err,error,*999)
                CASE(PROBLEM_SETUP_FINISH_ACTION)
                  !Finish the control loops
                  CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
                  CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,err,error,*999)
                  CALL CONTROL_LOOP_CREATE_FINISH(CONTROL_LOOP,err,error,*999)
                CASE DEFAULT
                  localError="The action type of "//TRIM(NumberToVString(PROBLEM_SETUP%ACTION_TYPE,"*",err,error))// &
                    & " for a setup type of "//TRIM(NumberToVString(PROBLEM_SETUP%SETUP_TYPE,"*",err,error))// &
                    & " is invalid for a ALE Stokes fluid."
                  CALL FlagError(localError,err,error,*999)
              END SELECT
            CASE(PROBLEM_SETUP_SOLVERS_TYPE)
              !Get the control loop
              CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
              CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,err,error,*999)
              SELECT CASE(PROBLEM_SETUP%ACTION_TYPE)
                CASE(PROBLEM_SETUP_START_ACTION)
                  !Start the solvers creation
                  CALL SOLVERS_CREATE_START(CONTROL_LOOP,SOLVERS,err,error,*999)
                  CALL SOLVERS_NUMBER_SET(SOLVERS,2,err,error,*999)
                  !Set the first solver to be a linear solver for the Laplace mesh movement problem
                  CALL SOLVERS_SOLVER_GET(SOLVERS,1,MESH_SOLVER,err,error,*999)
                  CALL SOLVER_TYPE_SET(MESH_SOLVER,SOLVER_LINEAR_TYPE,err,error,*999)
                  !Set solver defaults
                  CALL SOLVER_LIBRARY_TYPE_SET(MESH_SOLVER,SOLVER_PETSC_LIBRARY,err,error,*999)
                  !Set the solver to be a first order dynamic solver
                  CALL SOLVERS_SOLVER_GET(SOLVERS,2,SOLVER,err,error,*999)
                  CALL SOLVER_TYPE_SET(SOLVER,SOLVER_DYNAMIC_TYPE,err,error,*999)
                  CALL SOLVER_DYNAMIC_ORDER_SET(SOLVER,SOLVER_DYNAMIC_FIRST_ORDER,err,error,*999)
                  !Set solver defaults
                  CALL SOLVER_DYNAMIC_DEGREE_SET(SOLVER,SOLVER_DYNAMIC_FIRST_DEGREE,err,error,*999)
                  CALL SOLVER_DYNAMIC_SCHEME_SET(SOLVER,SOLVER_DYNAMIC_CRANK_NICOLSON_SCHEME,err,error,*999)
                  CALL SOLVER_LIBRARY_TYPE_SET(SOLVER,SOLVER_CMISS_LIBRARY,err,error,*999)
                CASE(PROBLEM_SETUP_FINISH_ACTION)
                  !Get the solvers
                  CALL CONTROL_LOOP_SOLVERS_GET(CONTROL_LOOP,SOLVERS,err,error,*999)
                  !Finish the solvers creation
                  CALL SOLVERS_CREATE_FINISH(SOLVERS,err,error,*999)
                CASE DEFAULT
                  localError="The action type of "//TRIM(NumberToVString(PROBLEM_SETUP%ACTION_TYPE,"*",err,error))// &
                   & " for a setup type of "//TRIM(NumberToVString(PROBLEM_SETUP%SETUP_TYPE,"*",err,error))// &
                   & " is invalid for a ALE Stokes fluid."
                  CALL FlagError(localError,err,error,*999)
              END SELECT
            CASE(PROBLEM_SETUP_SOLVER_EQUATIONS_TYPE)
              SELECT CASE(PROBLEM_SETUP%ACTION_TYPE)
                CASE(PROBLEM_SETUP_START_ACTION)
                  !Get the control loop
                  CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
                  CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,err,error,*999)
                  !Get the solver
                  CALL CONTROL_LOOP_SOLVERS_GET(CONTROL_LOOP,SOLVERS,err,error,*999)
                  CALL SOLVERS_SOLVER_GET(SOLVERS,1,MESH_SOLVER,err,error,*999)
                  !Create the solver equations
                  CALL SOLVER_EQUATIONS_CREATE_START(MESH_SOLVER,MESH_SOLVER_EQUATIONS,err,error,*999)
                  CALL SOLVER_EQUATIONS_LINEARITY_TYPE_SET(MESH_SOLVER_EQUATIONS,SOLVER_EQUATIONS_LINEAR,err,error,*999)
                  CALL SOLVER_EQUATIONS_TIME_DEPENDENCE_TYPE_SET(MESH_SOLVER_EQUATIONS,SOLVER_EQUATIONS_STATIC,err,error,*999)
                  CALL SOLVER_EQUATIONS_SPARSITY_TYPE_SET(MESH_SOLVER_EQUATIONS,SOLVER_SPARSE_MATRICES,err,error,*999)
                  CALL SOLVERS_SOLVER_GET(SOLVERS,2,SOLVER,err,error,*999)
                  !Create the solver equations
                  CALL SOLVER_EQUATIONS_CREATE_START(SOLVER,SOLVER_EQUATIONS,err,error,*999)
                  CALL SOLVER_EQUATIONS_LINEARITY_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_LINEAR,err,error,*999)
                  CALL SOLVER_EQUATIONS_TIME_DEPENDENCE_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_FIRST_ORDER_DYNAMIC,&
                    & err,error,*999)
                  CALL SOLVER_EQUATIONS_SPARSITY_TYPE_SET(SOLVER_EQUATIONS,SOLVER_SPARSE_MATRICES,err,error,*999)
                CASE(PROBLEM_SETUP_FINISH_ACTION)
                  !Get the control loop
                  CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
                  CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,err,error,*999)
                  !Get the solver equations
                  CALL CONTROL_LOOP_SOLVERS_GET(CONTROL_LOOP,SOLVERS,err,error,*999)
                  CALL SOLVERS_SOLVER_GET(SOLVERS,1,MESH_SOLVER,err,error,*999)
                  CALL SOLVER_SOLVER_EQUATIONS_GET(MESH_SOLVER,MESH_SOLVER_EQUATIONS,err,error,*999)
                  !Finish the solver equations creation
                  CALL SOLVER_EQUATIONS_CREATE_FINISH(MESH_SOLVER_EQUATIONS,err,error,*999)
                  CALL SOLVERS_SOLVER_GET(SOLVERS,2,SOLVER,err,error,*999)
                  CALL SOLVER_SOLVER_EQUATIONS_GET(SOLVER,SOLVER_EQUATIONS,err,error,*999)
                  !Finish the solver equations creation
                  CALL SOLVER_EQUATIONS_CREATE_FINISH(SOLVER_EQUATIONS,err,error,*999)
                CASE DEFAULT
                  localError="The action type of "//TRIM(NumberToVString(PROBLEM_SETUP%ACTION_TYPE,"*",err,error))// &
                    & " for a setup type of "//TRIM(NumberToVString(PROBLEM_SETUP%SETUP_TYPE,"*",err,error))// &
                    & " is invalid for a ALE Stokes fluid."
                  CALL FlagError(localError,err,error,*999)
              END SELECT
            CASE DEFAULT
              localError="The setup type of "//TRIM(NumberToVString(PROBLEM_SETUP%SETUP_TYPE,"*",err,error))// &
                & " is invalid for a ALE Stokes fluid."
              CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE DEFAULT
          localError="Problem subtype "//TRIM(NumberToVString(PROBLEM%SPECIFICATION(3),"*",err,error))// &
            & " is not valid for a Stokes fluid type of a fluid mechanics problem class."
          CALL FlagError(localError,err,error,*999)
      END SELECT
    ELSE
      CALL FlagError("Problem is not associated.",err,error,*999)
    ENDIF

    EXITS("STOKES_PROBLEM_SETUP")
    RETURN
999 ERRORSEXITS("STOKES_PROBLEM_SETUP",err,error)
    RETURN 1
  END SUBROUTINE STOKES_PROBLEM_SETUP

!
!================================================================================================================================
!

  !>Calculates the element stiffness matrices and RHS for a Stokes fluid finite element equations set.
  SUBROUTINE STOKES_FINITE_ELEMENT_CALCULATE(EQUATIONS_SET,ELEMENT_NUMBER,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to perform the finite element calculations on
    INTEGER(INTG), INTENT(IN) :: ELEMENT_NUMBER !<The element number to calculate
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) FIELD_VAR_TYPE,ng,mh,mhs,mi,ms,nh,nhs,ni,ns,MESH_COMPONENT1,MESH_COMPONENT2, nhs_max, mhs_max, nhs_min, mhs_min
    REAL(DP) :: JGW,SUM,DXI_DX(3,3),PHIMS,PHINS,MU_PARAM,RHO_PARAM,DPHIMS_DXI(3),DPHINS_DXI(3)
    LOGICAL :: updateStiffnessMatrix, updateDampingMatrix,updateRHSVector
    TYPE(BASIS_TYPE), POINTER :: DEPENDENT_BASIS,DEPENDENT_BASIS1,DEPENDENT_BASIS2,GEOMETRIC_BASIS,INDEPENDENT_BASIS
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(EquationsMappingLinearType), POINTER :: linearMapping
    TYPE(EquationsMappingDynamicType), POINTER :: dynamicMapping
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices
    TYPE(EquationsMatricesLinearType), POINTER :: linearMatrices
    TYPE(EquationsMatricesDynamicType), POINTER :: dynamicMatrices
    TYPE(EquationsMatricesRHSType), POINTER :: rhsVector
    TYPE(EquationsMatrixType), POINTER :: stiffnessMatrix, dampingMatrix
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(FIELD_TYPE), POINTER :: dependentField,geometricField,materialsField,independentField
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE
    TYPE(QUADRATURE_SCHEME_TYPE), POINTER :: QUADRATURE_SCHEME,QUADRATURE_SCHEME1,QUADRATURE_SCHEME2
    TYPE(VARYING_STRING) :: localError
    INTEGER:: xv,out
    REAL(DP) :: AG_MATRIX(256,256) ! "A" Matrix ("G"radient part) - maximum size allocated
    REAL(DP) :: AL_MATRIX(256,256) ! "A" Matrix ("L"aplace part) - maximum size allocated
    REAL(DP) :: BT_MATRIX(256,256) ! "B" "T"ranspose Matrix - maximum size allocated
    REAL(DP) :: MT_MATRIX(256,256) ! "M"ass "T"ime Matrix - maximum size allocated
    REAL(DP) :: CT_MATRIX(256,256) ! "C"onvective "T"erm Matrix - maximum size allocated
    REAL(DP) :: ALE_MATRIX(256,256) ! "A"rbitrary "L"agrangian "E"ulerian Matrix - maximum size allocated
    REAL(DP) :: RH_VECTOR(256) ! "R"ight "H"and vector - maximum size allocated
    REAL(DP) :: W_VALUE(3)
    REAL(DP)::  X(3)

!\todo: Reduce number of variables and parameters

    ENTERS("STOKES_FINITE_ELEMENT_CALCULATE",err,error,*999)

    out=0
    AG_MATRIX=0.0_DP
    AL_MATRIX=0.0_DP
    BT_MATRIX=0.0_DP
    MT_MATRIX=0.0_DP
    CT_MATRIX=0.0_DP
    ALE_MATRIX=0.0_DP
    RH_VECTOR=0.0_DP
    X=0.0_DP
!     L=10.0_DP

    updateStiffnessMatrix=.FALSE.
    updateDampingMatrix=.FALSE.
    updateRHSVector=.FALSE.

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(.NOT.ALLOCATED(EQUATIONS_SET%SPECIFICATION)) THEN
        CALL FlagError("Equations set specification is not allocated.",err,error,*999)
      ELSE IF(SIZE(EQUATIONS_SET%SPECIFICATION,1)/=3) THEN
        CALL FlagError("Equations set specification must have three entries for a Stokes type equations set.", &
          & err,error,*999)
      END IF
      EQUATIONS=>EQUATIONS_SET%EQUATIONS
      IF(ASSOCIATED(EQUATIONS)) THEN
        NULLIFY(vectorEquations)
        CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
        SELECT CASE(EQUATIONS_SET%SPECIFICATION(3))
          CASE(EQUATIONS_SET_STATIC_STOKES_SUBTYPE,EQUATIONS_SET_LAPLACE_STOKES_SUBTYPE, &
            & EQUATIONS_SET_TRANSIENT_STOKES_SUBTYPE,PROBLEM_ALE_STOKES_SUBTYPE,PROBLEM_PGM_STOKES_SUBTYPE)
            !Set Pointers
            dependentField=>equations%interpolation%dependentField
            independentField=>equations%interpolation%independentField
            geometricField=>equations%interpolation%geometricField
            materialsField=>equations%interpolation%materialsField
            vectorMatrices=>vectorEquations%vectorMatrices
            GEOMETRIC_BASIS=>geometricField%DECOMPOSITION%DOMAIN(geometricField%DECOMPOSITION%MESH_COMPONENT_NUMBER)%ptr% &
              & TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BASIS
            DEPENDENT_BASIS=>dependentField%DECOMPOSITION%DOMAIN(dependentField%DECOMPOSITION%MESH_COMPONENT_NUMBER)%ptr% &
              & TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BASIS
            QUADRATURE_SCHEME=>DEPENDENT_BASIS%QUADRATURE%QUADRATURE_SCHEME_MAP(BASIS_DEFAULT_QUADRATURE_SCHEME)%ptr
            rhsVector=>vectorMatrices%rhsVector
            vectorMapping=>vectorEquations%vectorMapping
            SELECT CASE(EQUATIONS_SET%SPECIFICATION(3))
              CASE(EQUATIONS_SET_STATIC_STOKES_SUBTYPE,EQUATIONS_SET_LAPLACE_STOKES_SUBTYPE)
                linearMatrices=>vectorMatrices%linearMatrices
                stiffnessMatrix=>linearMatrices%matrices(1)%ptr
                linearMapping=>vectorMapping%linearMapping
                FIELD_VARIABLE=>linearMapping%equationsMatrixToVarMaps(1)%variable
                stiffnessMatrix%elementMatrix%matrix=0.0_DP
                IF(ASSOCIATED(stiffnessMatrix)) updateStiffnessMatrix=stiffnessMatrix%updateMatrix
                IF(ASSOCIATED(rhsVector)) updateRHSVector=rhsVector%updateVector
              CASE(EQUATIONS_SET_TRANSIENT_STOKES_SUBTYPE)
                dynamicMatrices=>vectorMatrices%dynamicMatrices
                stiffnessMatrix=>dynamicMatrices%matrices(1)%ptr
                dampingMatrix=>dynamicMatrices%matrices(2)%ptr
                dynamicMapping=>vectorMapping%dynamicMapping
                FIELD_VARIABLE=>dynamicMapping%equationsMatrixToVarMaps(1)%variable
                FIELD_VAR_TYPE=FIELD_VARIABLE%VARIABLE_TYPE
                stiffnessMatrix%elementMatrix%matrix=0.0_DP
                dampingMatrix%elementMatrix%matrix=0.0_DP
                IF(ASSOCIATED(stiffnessMatrix)) updateStiffnessMatrix=stiffnessMatrix%updateMatrix
                IF(ASSOCIATED(dampingMatrix)) updateDampingMatrix=dampingMatrix%updateMatrix
                IF(ASSOCIATED(rhsVector)) updateRHSVector=rhsVector%updateVector
              CASE(EQUATIONS_SET_ALE_STOKES_SUBTYPE,EQUATIONS_SET_PGM_STOKES_SUBTYPE)
                independentField=>equations%interpolation%independentField
                INDEPENDENT_BASIS=>independentField%DECOMPOSITION%DOMAIN(independentField%DECOMPOSITION%MESH_COMPONENT_NUMBER)% &
                  & PTR%TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BASIS
                dynamicMatrices=>vectorMatrices%dynamicMatrices
                stiffnessMatrix=>dynamicMatrices%matrices(1)%ptr
                dampingMatrix=>dynamicMatrices%matrices(2)%ptr
                dynamicMapping=>vectorMapping%dynamicMapping
                FIELD_VARIABLE=>dynamicMapping%equationsMatrixToVarMaps(1)%variable
                FIELD_VAR_TYPE=FIELD_VARIABLE%VARIABLE_TYPE
                stiffnessMatrix%elementMatrix%matrix=0.0_DP
                dampingMatrix%elementMatrix%matrix=0.0_DP
                IF(ASSOCIATED(stiffnessMatrix)) updateStiffnessMatrix=stiffnessMatrix%updateMatrix
                IF(ASSOCIATED(dampingMatrix)) updateDampingMatrix=dampingMatrix%updateMatrix
                IF(ASSOCIATED(rhsVector)) updateRHSVector=rhsVector%updateVector
                CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_MESH_VELOCITY_SET_TYPE,ELEMENT_NUMBER, &
                  & equations%interpolation%independentInterpParameters(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
              CASE DEFAULT
                localError="Equations set subtype "//TRIM(NumberToVString(EQUATIONS_SET%SPECIFICATION(3),"*",err,error))// &
                  & " is not valid for a Stokes fluid type of a fluid mechanics equations set class."
                CALL FlagError(localError,err,error,*999)
            END SELECT
            CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER,equations%interpolation% &
              & geometricInterpParameters(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
            CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER,equations%interpolation% &
              & materialsInterpParameters(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)

            !Start looping over Gauss points
            DO ng=1,QUADRATURE_SCHEME%NUMBER_OF_GAUSS
              CALL FIELD_INTERPOLATE_GAUSS(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng,equations%interpolation% &
                & geometricInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
              CALL FIELD_INTERPOLATED_POINT_METRICS_CALCULATE(GEOMETRIC_BASIS%NUMBER_OF_XI,equations%interpolation% &
                & geometricInterpPointMetrics(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
              CALL FIELD_INTERPOLATE_GAUSS(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng,equations%interpolation% &
                & materialsInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
              IF(EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_ALE_STOKES_SUBTYPE.OR. &
                & EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_PGM_STOKES_SUBTYPE) THEN
                CALL FIELD_INTERPOLATE_GAUSS(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng,equations%interpolation% &
                  & independentInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
                  W_VALUE(1)=equations%interpolation%independentInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr%VALUES(1,NO_PART_DERIV)
                  W_VALUE(2)=equations%interpolation%independentInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr%VALUES(2,NO_PART_DERIV)
                  IF(FIELD_VARIABLE%NUMBER_OF_COMPONENTS==4) THEN
                    W_VALUE(3)=equations%interpolation%independentInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr%VALUES(3,NO_PART_DERIV)
                  END IF
              ELSE
                W_VALUE=0.0_DP
              END IF
              !Define MU_PARAM, viscosity=1
              MU_PARAM=equations%interpolation%materialsInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr%VALUES(1,NO_PART_DERIV)
              !Define RHO_PARAM, density=2
              RHO_PARAM=equations%interpolation%materialsInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr%VALUES(2,NO_PART_DERIV)
              !Calculate partial matrices
!\todo: Check time spent here
              IF(EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_STATIC_STOKES_SUBTYPE.OR. &
                & EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_LAPLACE_STOKES_SUBTYPE.OR. &
                & EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_ALE_STOKES_SUBTYPE.OR. &
                & EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_PGM_STOKES_SUBTYPE.OR. &
                & EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_TRANSIENT_STOKES_SUBTYPE) THEN
                !Loop over field components
                mhs=0
                DO mh=1,(FIELD_VARIABLE%NUMBER_OF_COMPONENTS-1)
                  MESH_COMPONENT1=FIELD_VARIABLE%COMPONENTS(mh)%MESH_COMPONENT_NUMBER
                  DEPENDENT_BASIS1=>dependentField%DECOMPOSITION%DOMAIN(MESH_COMPONENT1)%ptr% &
                    & TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BASIS
                  QUADRATURE_SCHEME1=>DEPENDENT_BASIS1%QUADRATURE%QUADRATURE_SCHEME_MAP(BASIS_DEFAULT_QUADRATURE_SCHEME)%ptr
                  JGW=equations%interpolation%geometricInterpPointMetrics(FIELD_U_VARIABLE_TYPE)%ptr%JACOBIAN* &
                    & QUADRATURE_SCHEME1%GAUSS_WEIGHTS(ng)
                  DO ms=1,DEPENDENT_BASIS1%NUMBER_OF_ELEMENT_PARAMETERS
                    mhs=mhs+1
                    nhs=0
                    IF(updateStiffnessMatrix.OR.updateDampingMatrix) THEN
                      !Loop over element columns
                      DO nh=1,(FIELD_VARIABLE%NUMBER_OF_COMPONENTS)

                        MESH_COMPONENT2=FIELD_VARIABLE%COMPONENTS(nh)%MESH_COMPONENT_NUMBER
                        DEPENDENT_BASIS2=>dependentField%DECOMPOSITION%DOMAIN(MESH_COMPONENT2)%ptr% &
                          & TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BASIS
                        QUADRATURE_SCHEME2=>DEPENDENT_BASIS2%QUADRATURE%QUADRATURE_SCHEME_MAP &
                          & (BASIS_DEFAULT_QUADRATURE_SCHEME)%ptr
                        ! JGW=equations%interpolation%geometricInterpPointMetrics%JACOBIAN*QUADRATURE_SCHEME2%&
                        ! &GAUSS_WEIGHTS(ng)
                        DO ns=1,DEPENDENT_BASIS2%NUMBER_OF_ELEMENT_PARAMETERS
                          nhs=nhs+1
                        !Calculate some variables used later on
                          DO ni=1,DEPENDENT_BASIS2%NUMBER_OF_XI
                            DO mi=1,DEPENDENT_BASIS1%NUMBER_OF_XI
                              DXI_DX(mi,ni)=equations%interpolation%geometricInterpPointMetrics(FIELD_U_VARIABLE_TYPE)%ptr% &
                                & DXI_DX(mi,ni)
                            END DO
                            DPHIMS_DXI(ni)=QUADRATURE_SCHEME1%GAUSS_BASIS_FNS(ms,PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(ni),ng)
                            DPHINS_DXI(ni)=QUADRATURE_SCHEME2%GAUSS_BASIS_FNS(ns,PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(ni),ng)
                          END DO !ni
                          PHIMS=QUADRATURE_SCHEME1%GAUSS_BASIS_FNS(ms,NO_PART_DERIV,ng)
                          PHINS=QUADRATURE_SCHEME2%GAUSS_BASIS_FNS(ns,NO_PART_DERIV,ng)
                        !                         DO mi=1,DEPENDENT_BASIS1%NUMBER_OF_XI
                        !                           DO ni=1,DEPENDENT_BASIS2%NUMBER_OF_XI
                        !                             SUM=SUM-MU_PARAM*DPHIMSS_DXI(mi)*DPHINSS_DXI(ni)*equations%interpolation%geometricInterpPointMetrics%GU(mi,ni)
                        !                           ENDDO !ni
                        !                         ENDDO !mi
                          IF(updateStiffnessMatrix) THEN

                            !LAPLACE TYPE
                            IF(nh==mh) THEN
                              SUM=0.0_DP
                              !Calculate SUM
                              DO xv=1,DEPENDENT_BASIS1%NUMBER_OF_XI
                                DO mi=1,DEPENDENT_BASIS1%NUMBER_OF_XI
                                  DO ni=1,DEPENDENT_BASIS2%NUMBER_OF_XI
                                    SUM=SUM+MU_PARAM*DPHINS_DXI(ni)*DXI_DX(ni,xv)*DPHIMS_DXI(mi)*DXI_DX(mi,xv)
                                  ENDDO !ni
                                ENDDO !mi
                              ENDDO !x
                              !Calculate MATRIX
                              AL_MATRIX(mhs,nhs)=AL_MATRIX(mhs,nhs)+SUM*JGW
                            END IF

                          END IF
                          !Calculate standard matrix (gradient transpose type)
                          IF(updateStiffnessMatrix) THEN

                            IF(EQUATIONS_SET%SPECIFICATION(3)/=EQUATIONS_SET_LAPLACE_STOKES_SUBTYPE) THEN
                              IF(nh<FIELD_VARIABLE%NUMBER_OF_COMPONENTS) THEN
                                SUM=0.0_DP
                                !Calculate SUM
                                DO mi=1,DEPENDENT_BASIS1%NUMBER_OF_XI
                                  DO ni=1,DEPENDENT_BASIS2%NUMBER_OF_XI
                                    !note mh/nh derivative in DXI_DX
                                    SUM=SUM+MU_PARAM*DPHINS_DXI(mi)*DXI_DX(mi,mh)*DPHIMS_DXI(ni)*DXI_DX(ni,nh)
                                  ENDDO !ni
                                ENDDO !mi
                                !Calculate MATRIX
                                AG_MATRIX(mhs,nhs)=AG_MATRIX(mhs,nhs)+SUM*JGW
                              END IF
                            END IF

                          END IF
                          !Calculate ALE matric contribution
                          IF(updateStiffnessMatrix) THEN

                            IF(EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_ALE_STOKES_SUBTYPE.OR. &
                              & EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_PGM_STOKES_SUBTYPE) THEN
                              IF(nh==mh) THEN
                                SUM=0.0_DP
                                !Calculate SUM
                                DO mi=1,DEPENDENT_BASIS1%NUMBER_OF_XI
                                  DO ni=1,DEPENDENT_BASIS1%NUMBER_OF_XI
                                    SUM=SUM-RHO_PARAM*W_VALUE(mi)*DPHINS_DXI(ni)*DXI_DX(ni,mi)*PHIMS
                                  ENDDO !ni
                                ENDDO !mi
                                !Calculate MATRIX
                                ALE_MATRIX(mhs,nhs)=ALE_MATRIX(mhs,nhs)+SUM*JGW
                              END IF
                            END IF

                          END IF
                          !Calculate pressure contribution (B transpose type)
                          IF(updateStiffnessMatrix) THEN

                            !LAPLACE TYPE
                            IF(nh==FIELD_VARIABLE%NUMBER_OF_COMPONENTS) THEN
                              SUM=0.0_DP
                              !Calculate SUM
                              DO ni=1,DEPENDENT_BASIS1%NUMBER_OF_XI
                                SUM=SUM-PHINS*DPHIMS_DXI(ni)*DXI_DX(ni,mh)
                              ENDDO !ni
                              !Calculate MATRIX
                              BT_MATRIX(mhs,nhs)=BT_MATRIX(mhs,nhs)+SUM*JGW
                            END IF

                          END IF
                          !Calculate mass matrix if needed
                          IF(EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_TRANSIENT_STOKES_SUBTYPE.OR. &
                            & EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_ALE_STOKES_SUBTYPE.OR. &
                            & EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_PGM_STOKES_SUBTYPE) THEN
                            IF(updateDampingMatrix) THEN
                              IF(nh==mh) THEN
                                SUM=0.0_DP
                                !Calculate SUM
                                SUM=PHIMS*PHINS*RHO_PARAM
                                !Calculate MATRIX
                                MT_MATRIX(mhs,nhs)=MT_MATRIX(mhs,nhs)+SUM*JGW
                              END IF
                            END IF
                          END IF
                        ENDDO !ns
                      ENDDO !nh
                    ENDIF
                  ENDDO !ms
                ENDDO !mh

                !Calculate analytic RHS
                IF(ASSOCIATED(EQUATIONS_SET%ANALYTIC)) THEN
                  IF(EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE==EQUATIONS_SET_STOKES_EQUATION_TWO_DIM_1.OR. &
                    & EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE==EQUATIONS_SET_STOKES_EQUATION_TWO_DIM_2.OR. &
                    & EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE==EQUATIONS_SET_STOKES_EQUATION_TWO_DIM_3.OR. &
                    & EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE==EQUATIONS_SET_STOKES_EQUATION_TWO_DIM_4.OR. &
                    & EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE==EQUATIONS_SET_STOKES_EQUATION_TWO_DIM_5.OR. &
                    & EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE==EQUATIONS_SET_STOKES_EQUATION_THREE_DIM_1.OR. &
                    & EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE==EQUATIONS_SET_STOKES_EQUATION_THREE_DIM_2.OR. &
                    & EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE==EQUATIONS_SET_STOKES_EQUATION_THREE_DIM_3.OR. &
                    & EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE==EQUATIONS_SET_STOKES_EQUATION_THREE_DIM_4.OR. &
                    & EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE==EQUATIONS_SET_STOKES_EQUATION_THREE_DIM_5) THEN

                    mhs=0
                    DO mh=1,(FIELD_VARIABLE%NUMBER_OF_COMPONENTS-1)
                      MESH_COMPONENT1=FIELD_VARIABLE%COMPONENTS(mh)%MESH_COMPONENT_NUMBER
                      DEPENDENT_BASIS1=>dependentField%DECOMPOSITION%DOMAIN(MESH_COMPONENT1)%ptr% &
                        & TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BASIS
                      QUADRATURE_SCHEME1=>DEPENDENT_BASIS1%QUADRATURE%QUADRATURE_SCHEME_MAP(BASIS_DEFAULT_QUADRATURE_SCHEME)%ptr
                      JGW=equations%interpolation%geometricInterpPointMetrics(FIELD_U_VARIABLE_TYPE)%ptr%JACOBIAN* &
                        & QUADRATURE_SCHEME1%GAUSS_WEIGHTS(ng)
                      DO ms=1,DEPENDENT_BASIS1%NUMBER_OF_ELEMENT_PARAMETERS
                        mhs=mhs+1
                        PHIMS=QUADRATURE_SCHEME1%GAUSS_BASIS_FNS(ms,NO_PART_DERIV,ng)
                        !note mh value derivative
                        SUM=0.0_DP
                        X(1) = equations%interpolation%geometricInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr%VALUES(1,1)
                        X(2) = equations%interpolation%geometricInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr%VALUES(2,1)
                        IF(DEPENDENT_BASIS1%NUMBER_OF_XI==3) THEN
                          X(3) = equations%interpolation%geometricInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr%VALUES(3,1)
                        END IF
                        IF(EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE==EQUATIONS_SET_STOKES_EQUATION_TWO_DIM_1) THEN
                          IF(mh==1) THEN
                            !Calculate SUM
                            SUM=0.0_DP
                          ELSE IF(mh==2) THEN
                            !Calculate SUM
                            SUM=PHIMS*(-2.0_DP*MU_PARAM/10.0_DP**2)
                          ENDIF
                        ELSE IF(EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE==EQUATIONS_SET_STOKES_EQUATION_TWO_DIM_2) THEN
                          IF(mh==1) THEN
                            !Calculate SUM
                            SUM=0.0_DP
                          ELSE IF(mh==2) THEN
                            !Calculate SUM
                            SUM=PHIMS*(-4.0_DP*MU_PARAM/100.0_DP*EXP((X(1)-X(2))/10.0_DP))
                          ENDIF
                        ELSE IF(EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE==EQUATIONS_SET_STOKES_EQUATION_TWO_DIM_3) THEN
                          IF(mh==1) THEN
                            !Calculate SUM
                            SUM=0.0_DP
                          ELSE IF(mh==2) THEN
                            !Calculate SUM
                            SUM=PHIMS*(16.0_DP*MU_PARAM*PI*PI/100.0_DP*COS(2.0_DP*PI*X(2)/10.0_DP)*COS(2.0_DP*PI*X(1)/10.0_DP))
                          ENDIF
                        ELSE IF(EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE==EQUATIONS_SET_STOKES_EQUATION_TWO_DIM_4) THEN
!                           do nothing!
                        ELSE IF(EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE==EQUATIONS_SET_STOKES_EQUATION_TWO_DIM_5) THEN
!                           do nothing!
                        ELSE IF(EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE==EQUATIONS_SET_STOKES_EQUATION_THREE_DIM_4) THEN
!                           do nothing!
                        ELSE IF(EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE==EQUATIONS_SET_STOKES_EQUATION_THREE_DIM_5) THEN
!                           do nothing!
                        ELSE IF(EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE==EQUATIONS_SET_STOKES_EQUATION_THREE_DIM_1) THEN
                          IF(mh==1) THEN
                            !Calculate SUM
                            SUM=0.0_DP
                          ELSE IF(mh==2) THEN
                            !Calculate SUM
                            SUM=PHIMS*(-4.0_DP*MU_PARAM/100.0_DP)
                          ELSE IF(mh==3) THEN
                            !Calculate SUM
                            SUM=PHIMS*(-4.0_DP*MU_PARAM/100.0_DP)
                          ENDIF
                        ELSE IF(EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE==EQUATIONS_SET_STOKES_EQUATION_THREE_DIM_2) THEN
                          IF(mh==1) THEN
                            !Calculate SUM
                            SUM=0.0_DP
                          ELSE IF(mh==2) THEN
                            !Calculate SUM
                            SUM=PHIMS*(-2.0_DP*MU_PARAM/100.0_DP*(2.0_DP*EXP((X(1)-X(2))/10.0_DP)+EXP((X(2)-X(3))/10.0_DP)))
                          ELSE IF(mh==3) THEN
                            !Calculate SUM
                            SUM=PHIMS*(-2.0_DP*MU_PARAM/100.0_DP*(2.0_DP*EXP((X(3)-X(1))/10.0_DP)+EXP((X(2)-X(3))/10.0_DP)))
                          ENDIF
                        ELSE IF(EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE==EQUATIONS_SET_STOKES_EQUATION_THREE_DIM_3) THEN
                          IF(mh==1) THEN
                            !Calculate SUM
                            SUM=0.0_DP
                          ELSE IF(mh==2) THEN
                            !Calculate SUM
                            SUM=PHIMS*(36*MU_PARAM*PI**2/100.0_DP*COS(2.0_DP*PI*X(2)/10.0_DP)*SIN(2.0_DP*PI*X(3)/10.0_DP)* &
                              & COS(2.0_DP*PI*X(1)/10.0_DP))
                          ELSE IF(mh==3) THEN
                            !Calculate SUM
                            SUM=0.0_DP
                          ENDIF
                        ENDIF
                        !Calculate RH VECTOR
                        RH_VECTOR(mhs)=RH_VECTOR(mhs)+SUM*JGW
                      ENDDO !ms
                    ENDDO !mh
                  ELSE
                    RH_VECTOR(mhs)=0.0_DP
                  ENDIF
                ENDIF
              END IF
            ENDDO !ng
            !Assemble matrices calculated above
            mhs_min=mhs
            mhs_max=nhs
            nhs_min=mhs
            nhs_max=nhs
            IF(EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_STATIC_STOKES_SUBTYPE.OR.  &
              & EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_LAPLACE_STOKES_SUBTYPE.OR. &
              & EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_ALE_STOKES_SUBTYPE.OR. &
              & EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_PGM_STOKES_SUBTYPE.OR. &
              & EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_TRANSIENT_STOKES_SUBTYPE) THEN
              IF(updateStiffnessMatrix) THEN
                stiffnessMatrix%elementMatrix%matrix(1:mhs_min,1:nhs_min)=AL_MATRIX(1:mhs_min,1:nhs_min)+AG_MATRIX(1:mhs_min, &
                  & 1:nhs_min)+ALE_MATRIX(1:mhs_min,1:nhs_min)
                stiffnessMatrix%elementMatrix%matrix(1:mhs_min,nhs_min+1:nhs_max)=BT_MATRIX(1:mhs_min,nhs_min+1:nhs_max)
                DO mhs=mhs_min+1,mhs_max
                  DO nhs=1,nhs_min
                    !Transpose pressure type entries for mass equation
                    stiffnessMatrix%elementMatrix%matrix(mhs,nhs)=stiffnessMatrix%elementMatrix%matrix(nhs,mhs)
                  END DO
                END DO
              END IF
            END IF
            IF(EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_TRANSIENT_STOKES_SUBTYPE.OR. &
              & EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_ALE_STOKES_SUBTYPE.OR. &
              & EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_PGM_STOKES_SUBTYPE) THEN
              IF(updateDampingMatrix) THEN
                dampingMatrix%elementMatrix%matrix(1:mhs_min,1:nhs_min)=MT_MATRIX(1:mhs_min,1:nhs_min)
              END IF
            END IF
          !Assemble RHS vector
          IF(rhsVector%firstAssembly) THEN
            IF(updateRHSVector) THEN
              rhsVector%elementVector%vector(1:mhs_max)=RH_VECTOR(1:mhs_max)
            ENDIF
          ENDIF
          !Scale factor adjustment
            IF(dependentField%SCALINGS%SCALING_TYPE/=FIELD_NO_SCALING) THEN
              CALL Field_InterpolationParametersScaleFactorsElementGet(ELEMENT_NUMBER,equations%interpolation% &
                & dependentInterpParameters(FIELD_VAR_TYPE)%ptr,err,error,*999)
              mhs=0
              DO mh=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
                !Loop over element rows
                MESH_COMPONENT1=FIELD_VARIABLE%COMPONENTS(mh)%MESH_COMPONENT_NUMBER
                DEPENDENT_BASIS1=>dependentField%DECOMPOSITION%DOMAIN(MESH_COMPONENT1)%ptr% &
                  & TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BASIS
                DO ms=1,DEPENDENT_BASIS1%NUMBER_OF_ELEMENT_PARAMETERS
                  mhs=mhs+1
                  nhs=0
                   IF(updateStiffnessMatrix.OR.updateDampingMatrix) THEN
                    !Loop over element columns
                    DO nh=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
                      MESH_COMPONENT2=FIELD_VARIABLE%COMPONENTS(nh)%MESH_COMPONENT_NUMBER
                      DEPENDENT_BASIS2=>dependentField%DECOMPOSITION%DOMAIN(MESH_COMPONENT2)%ptr% &
                        & TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BASIS
                      DO ns=1,DEPENDENT_BASIS2%NUMBER_OF_ELEMENT_PARAMETERS
                        nhs=nhs+1
                        IF(updateStiffnessMatrix)THEN
                          stiffnessMatrix%elementMatrix%matrix(mhs,nhs)=stiffnessMatrix%elementMatrix%matrix(mhs,nhs)* &
                            & equations%interpolation%dependentInterpParameters(FIELD_VAR_TYPE)%ptr%SCALE_FACTORS(ms,mh)* &
                            & equations%interpolation%dependentInterpParameters(FIELD_VAR_TYPE)%ptr%SCALE_FACTORS(ns,nh)
                        END IF
                        IF(updateDampingMatrix)THEN
                          dampingMatrix%elementMatrix%matrix(mhs,nhs)=dampingMatrix%elementMatrix%matrix(mhs,nhs)* &
                            & equations%interpolation%dependentInterpParameters(FIELD_VAR_TYPE)%ptr%SCALE_FACTORS(ms,mh)* &
                            & equations%interpolation%dependentInterpParameters(FIELD_VAR_TYPE)%ptr%SCALE_FACTORS(ns,nh)
                        END IF
                      ENDDO !ns
                    ENDDO !nh
                  ENDIF
                  IF(updateRHSVector) rhsVector%elementVector%vector(mhs)=rhsVector%elementVector%vector(mhs)* &
                    & equations%interpolation%dependentInterpParameters(FIELD_VAR_TYPE)%ptr%SCALE_FACTORS(ms,mh)
                ENDDO !ms
              ENDDO !mh
            ENDIF
          CASE DEFAULT
            localError="Equations set subtype "//TRIM(NumberToVString(EQUATIONS_SET%SPECIFICATION(3),"*",err,error))// &
              & " is not valid for a Stokes fluid type of a fluid mechanics equations set class."
            CALL FlagError(localError,err,error,*999)
        END SELECT
      ELSE
        CALL FlagError("Equations set equations is not associated.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Equations set is not associated.",err,error,*999)
    ENDIF

    EXITS("STOKES_FINITE_ELEMENT_CALCULATE")
    RETURN
999 ERRORSEXITS("STOKES_FINITE_ELEMENT_CALCULATE",err,error)
    RETURN 1
  END SUBROUTINE STOKES_FINITE_ELEMENT_CALCULATE

  !
  !================================================================================================================================
  !

  !>Sets up the Stokes problem post solve.
  SUBROUTINE STOKES_POST_SOLVE(CONTROL_LOOP,SOLVER,err,error,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP !<A pointer to the control loop to solve.
    TYPE(SOLVER_TYPE), POINTER :: SOLVER!<A pointer to the solver
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(SOLVER_TYPE), POINTER :: SOLVER2 !<A pointer to the solver
    TYPE(VARYING_STRING) :: localError

    ENTERS("STOKES_POST_SOLVE",err,error,*999)
    NULLIFY(SOLVER2)

    IF(ASSOCIATED(CONTROL_LOOP)) THEN
      IF(ASSOCIATED(SOLVER)) THEN
        IF(ASSOCIATED(CONTROL_LOOP%PROBLEM)) THEN
          IF(.NOT.ALLOCATED(control_loop%problem%specification)) THEN
            CALL FlagError("Problem specification is not allocated.",err,error,*999)
          ELSE IF(SIZE(control_loop%problem%specification,1)<3) THEN
            CALL FlagError("Problem specification must have three entries for a Stokes problem.",err,error,*999)
          END IF
          SELECT CASE(CONTROL_LOOP%PROBLEM%SPECIFICATION(3))
            CASE(PROBLEM_STATIC_STOKES_SUBTYPE,PROBLEM_LAPLACE_STOKES_SUBTYPE)
              CALL STOKES_POST_SOLVE_OUTPUT_DATA(CONTROL_LOOP,SOLVER,err,error,*999)
            CASE(PROBLEM_PGM_STOKES_SUBTYPE)
              CALL STOKES_POST_SOLVE_OUTPUT_DATA(CONTROL_LOOP,SOLVER,err,error,*999)
            CASE(PROBLEM_TRANSIENT_STOKES_SUBTYPE)
              CALL STOKES_POST_SOLVE_OUTPUT_DATA(CONTROL_LOOP,SOLVER,err,error,*999)
            CASE(PROBLEM_ALE_STOKES_SUBTYPE)
              !Post solve for the linear solver
              IF(SOLVER%SOLVE_TYPE==SOLVER_LINEAR_TYPE) THEN
                CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"Mesh movement post solve... ",err,error,*999)
                CALL SOLVERS_SOLVER_GET(SOLVER%SOLVERS,2,SOLVER2,err,error,*999)
                IF(ASSOCIATED(SOLVER2%DYNAMIC_SOLVER)) THEN
                  SOLVER2%DYNAMIC_SOLVER%ALE=.TRUE.
                ELSE
                  CALL FlagError("Dynamic solver is not associated for ALE problem.",err,error,*999)
                END IF
              !Post solve for the linear solver
              ELSE IF(SOLVER%SOLVE_TYPE==SOLVER_DYNAMIC_TYPE) THEN
                CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"ALE Stokes post solve... ",err,error,*999)
                CALL STOKES_POST_SOLVE_OUTPUT_DATA(CONTROL_LOOP,SOLVER,err,error,*999)
              END IF
            CASE DEFAULT
              localError="Problem subtype "//TRIM(NumberToVString(CONTROL_LOOP%PROBLEM%SPECIFICATION(3),"*",err,error))// &
                & " is not valid for a Stokes fluid type of a fluid mechanics problem class."
              CALL FlagError(localError,err,error,*999)
          END SELECT
        ELSE
          CALL FlagError("Problem is not associated.",err,error,*999)
        ENDIF
      ELSE
        CALL FlagError("Solver is not associated.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Control loop is not associated.",err,error,*999)
    ENDIF

    EXITS("STOKES_POST_SOLVE")
    RETURN
999 ERRORSEXITS("STOKES_POST_SOLVE",err,error)
    RETURN 1
  END SUBROUTINE STOKES_POST_SOLVE

  !
  !================================================================================================================================
  !

  !>Sets up the Stokes problem pre solve.
  SUBROUTINE STOKES_PRE_SOLVE(CONTROL_LOOP,SOLVER,err,error,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP !<A pointer to the control loop to solve.
    TYPE(SOLVER_TYPE), POINTER :: SOLVER !<A pointer to the solver
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(SOLVER_TYPE), POINTER :: SOLVER2 !<A pointer to the solvers
    TYPE(VARYING_STRING) :: localError

    ENTERS("STOKES_PRE_SOLVE",err,error,*999)
    NULLIFY(SOLVER2)

    IF(ASSOCIATED(CONTROL_LOOP)) THEN
      IF(ASSOCIATED(SOLVER)) THEN
        IF(ASSOCIATED(CONTROL_LOOP%PROBLEM)) THEN
          IF(.NOT.ALLOCATED(control_loop%problem%specification)) THEN
            CALL FlagError("Problem specification is not allocated.",err,error,*999)
          ELSE IF(SIZE(control_loop%problem%specification,1)<3) THEN
            CALL FlagError("Problem specification must have three entries for a Stokes problem.",err,error,*999)
          END IF
          SELECT CASE(CONTROL_LOOP%PROBLEM%SPECIFICATION(3))
            CASE(PROBLEM_STATIC_STOKES_SUBTYPE,PROBLEM_LAPLACE_STOKES_SUBTYPE)
              ! do nothing ???
            CASE(PROBLEM_TRANSIENT_STOKES_SUBTYPE)
              ! do nothing ???
                CALL STOKES_PRE_SOLVE_UPDATE_BOUNDARY_CONDITIONS(CONTROL_LOOP,SOLVER,err,error,*999)
            CASE(PROBLEM_PGM_STOKES_SUBTYPE)
              ! do nothing ???
              !First update mesh and calculates boundary velocity values
              CALL STOKES_PRE_SOLVE_ALE_UPDATE_MESH(CONTROL_LOOP,SOLVER,err,error,*999)
              !Then apply both normal and moving mesh boundary conditions
              CALL STOKES_PRE_SOLVE_UPDATE_BOUNDARY_CONDITIONS(CONTROL_LOOP,SOLVER,err,error,*999)
            CASE(PROBLEM_ALE_STOKES_SUBTYPE)
              !Pre solve for the linear solver
              IF(SOLVER%SOLVE_TYPE==SOLVER_LINEAR_TYPE) THEN
                CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"Mesh movement pre solve... ",err,error,*999)
                !Update boundary conditions for mesh-movement
                CALL STOKES_PRE_SOLVE_UPDATE_BOUNDARY_CONDITIONS(CONTROL_LOOP,SOLVER,err,error,*999)
                CALL SOLVERS_SOLVER_GET(SOLVER%SOLVERS,2,SOLVER2,err,error,*999)
                IF(ASSOCIATED(SOLVER2%DYNAMIC_SOLVER)) THEN
!\todo: Avoid ALE flag in future
                  SOLVER2%DYNAMIC_SOLVER%ALE=.FALSE.
                ELSE
                  CALL FlagError("Dynamic solver is not associated for ALE problem.",err,error,*999)
                END IF
                !Update material properties for Laplace mesh movement
                CALL STOKES_PRE_SOLVE_ALE_UPDATE_PARAMETERS(CONTROL_LOOP,SOLVER,err,error,*999)
              !Pre solve for the linear solver
              ELSE IF(SOLVER%SOLVE_TYPE==SOLVER_DYNAMIC_TYPE) THEN
                CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"ALE Stokes pre solve... ",err,error,*999)
                IF(SOLVER%DYNAMIC_SOLVER%ALE) THEN
                  !First update mesh and calculates boundary velocity values
                  CALL STOKES_PRE_SOLVE_ALE_UPDATE_MESH(CONTROL_LOOP,SOLVER,err,error,*999)
                  !Then apply both normal and moving mesh boundary conditions
                  CALL STOKES_PRE_SOLVE_UPDATE_BOUNDARY_CONDITIONS(CONTROL_LOOP,SOLVER,err,error,*999)
                ELSE
                  CALL FlagError("Mesh motion calculation not successful for ALE problem.",err,error,*999)
                END IF
              ELSE
                CALL FlagError("Solver type is not associated for ALE problem.",err,error,*999)
              END IF
            CASE DEFAULT
              localError="Problem subtype "//TRIM(NumberToVString(CONTROL_LOOP%PROBLEM%SPECIFICATION(3),"*",err,error))// &
                & " is not valid for a Stokes fluid type of a fluid mechanics problem class."
              CALL FlagError(localError,err,error,*999)
          END SELECT
        ELSE
          CALL FlagError("Problem is not associated.",err,error,*999)
        ENDIF
      ELSE
        CALL FlagError("Solver is not associated.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Control loop is not associated.",err,error,*999)
    ENDIF

    EXITS("STOKES_PRE_SOLVE")
    RETURN
999 ERRORSEXITS("STOKES_PRE_SOLVE",err,error)
    RETURN 1
  END SUBROUTINE STOKES_PRE_SOLVE
  !
  !================================================================================================================================
  !

  !>Update boundary conditions for Stokes flow pre solve
  SUBROUTINE STOKES_PRE_SOLVE_UPDATE_BOUNDARY_CONDITIONS(CONTROL_LOOP,SOLVER,err,error,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP !<A pointer to the control loop to solve.
    TYPE(SOLVER_TYPE), POINTER :: SOLVER !<A pointer to the solver
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS  !<A pointer to the solver equations
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: SOLVER_MAPPING !<A pointer to the solver mapping
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set
    TYPE(EquationsType), POINTER :: EQUATIONS
    TYPE(VARYING_STRING) :: localError
    TYPE(BOUNDARY_CONDITIONS_VARIABLE_TYPE), POINTER :: BOUNDARY_CONDITIONS_VARIABLE
    TYPE(BOUNDARY_CONDITIONS_TYPE), POINTER :: BOUNDARY_CONDITIONS
    TYPE(FIELD_TYPE), POINTER :: dependentField,geometricField,materialsField
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE,GEOMETRIC_VARIABLE
    TYPE(DOMAIN_TYPE), POINTER :: DOMAIN
    TYPE(DOMAIN_NODES_TYPE), POINTER :: DOMAIN_NODES
    TYPE(FIELD_INTERPOLATED_POINT_PTR_TYPE), POINTER :: INTERPOLATED_POINT(:)
    TYPE(FIELD_INTERPOLATION_PARAMETERS_PTR_TYPE), POINTER :: INTERPOLATION_PARAMETERS(:)
    REAL(DP) :: CURRENT_TIME,TIME_INCREMENT,DISPLACEMENT_VALUE,VALUE,XI_COORDINATES(3)
    REAL(DP) :: T_COORDINATES(20,3)
    INTEGER(INTG) :: NUMBER_OF_DIMENSIONS,BOUNDARY_CONDITION_CHECK_VARIABLE,GLOBAL_DERIV_INDEX,node_idx,variable_type
    INTEGER(INTG) :: variable_idx,local_ny,ANALYTIC_FUNCTION_TYPE,component_idx,deriv_idx,dim_idx
    INTEGER(INTG) :: element_idx,en_idx,I,J,K,number_of_nodes_xic(3)
    REAL(DP) :: X(3),MU_PARAM,RHO_PARAM
    REAL(DP), POINTER :: MESH_VELOCITY_VALUES(:), GEOMETRIC_PARAMETERS(:)
    REAL(DP), POINTER :: BOUNDARY_VALUES(:)


    ENTERS("STOKES_PRE_SOLVE_UPDATE_BOUNDARY_CONDITIONS",err,error,*999)

    IF(ASSOCIATED(CONTROL_LOOP)) THEN
      CALL CONTROL_LOOP_CURRENT_TIMES_GET(CONTROL_LOOP,CURRENT_TIME,TIME_INCREMENT,err,error,*999)
!       write(*,*)'CURRENT_TIME = ',CURRENT_TIME
!       write(*,*)'TIME_INCREMENT = ',TIME_INCREMENT
      IF(ASSOCIATED(SOLVER)) THEN
        IF(ASSOCIATED(CONTROL_LOOP%PROBLEM)) THEN
          IF(.NOT.ALLOCATED(control_loop%problem%specification)) THEN
            CALL FlagError("Problem specification is not allocated.",err,error,*999)
          ELSE IF(SIZE(control_loop%problem%specification,1)<3) THEN
            CALL FlagError("Problem specification must have three entries for a Stokes problem.",err,error,*999)
          END IF
          SELECT CASE(CONTROL_LOOP%PROBLEM%SPECIFICATION(3))
            CASE(PROBLEM_STATIC_STOKES_SUBTYPE,PROBLEM_LAPLACE_STOKES_SUBTYPE)
              ! do nothing ???
            CASE(PROBLEM_TRANSIENT_STOKES_SUBTYPE)
              SOLVER_EQUATIONS=>SOLVER%SOLVER_EQUATIONS
              IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
                SOLVER_MAPPING=>SOLVER_equations%SOLVER_MAPPING
                EQUATIONS=>SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(1)%EQUATIONS
                IF(ASSOCIATED(EQUATIONS)) THEN
                  EQUATIONS_SET=>equations%equationsSet
                  IF(ASSOCIATED(EQUATIONS_SET%ANALYTIC)) THEN
                    IF(EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE==EQUATIONS_SET_STOKES_EQUATION_TWO_DIM_4 .OR. &
                      & EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE==EQUATIONS_SET_STOKES_EQUATION_THREE_DIM_1 .OR. &
                      & EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE==EQUATIONS_SET_STOKES_EQUATION_THREE_DIM_4 .OR. &
                      & EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE==EQUATIONS_SET_STOKES_EQUATION_THREE_DIM_5 .OR. &
                      & EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE==EQUATIONS_SET_STOKES_EQUATION_TWO_DIM_5) THEN
                      IF(ASSOCIATED(EQUATIONS_SET)) THEN
                        IF(ASSOCIATED(EQUATIONS_SET%ANALYTIC)) THEN
                          dependentField=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD
                          IF(ASSOCIATED(dependentField)) THEN
                            geometricField=>EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD
                            IF(ASSOCIATED(geometricField)) THEN
                              CALL FIELD_NUMBER_OF_COMPONENTS_GET(geometricField,FIELD_U_VARIABLE_TYPE,&
                                & NUMBER_OF_DIMENSIONS,err,error,*999)
                              NULLIFY(INTERPOLATION_PARAMETERS)
                              NULLIFY(INTERPOLATED_POINT)
                              CALL FIELD_INTERPOLATION_PARAMETERS_INITIALISE(geometricField,INTERPOLATION_PARAMETERS,err,error, &
                                & *999)
                              CALL FIELD_INTERPOLATED_POINTS_INITIALISE(INTERPOLATION_PARAMETERS,INTERPOLATED_POINT,err,error,*999)
                              NULLIFY(GEOMETRIC_VARIABLE)
                              CALL Field_VariableGet(geometricField,FIELD_U_VARIABLE_TYPE,GEOMETRIC_VARIABLE,err,error,*999)
                              NULLIFY(GEOMETRIC_PARAMETERS)
                              CALL FIELD_PARAMETER_SET_DATA_GET(geometricField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,&
                                & GEOMETRIC_PARAMETERS,err,error,*999)
                               DO variable_idx=1,dependentField%NUMBER_OF_VARIABLES
                                variable_type=dependentField%VARIABLES(variable_idx)%variable_TYPE
                                FIELD_VARIABLE=>dependentField%VARIABLE_TYPE_MAP(variable_type)%ptr
                                IF(ASSOCIATED(FIELD_VARIABLE)) THEN
                                  DO component_idx=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
                                    IF(FIELD_VARIABLE%COMPONENTS(component_idx)%INTERPOLATION_TYPE== &
                                      & FIELD_NODE_BASED_INTERPOLATION) THEN
                                      DOMAIN=>FIELD_VARIABLE%COMPONENTS(component_idx)%DOMAIN
                                      IF(ASSOCIATED(DOMAIN)) THEN
                                        IF(ASSOCIATED(DOMAIN%TOPOLOGY)) THEN
                                          DOMAIN_NODES=>DOMAIN%TOPOLOGY%NODES
                                          IF(ASSOCIATED(DOMAIN_NODES)) THEN
                                            !Should be replaced by boundary node flag
                                            DO node_idx=1,DOMAIN_NODES%NUMBER_OF_NODES
                                              element_idx=DOMAIN%topology%nodes%nodes(node_idx)%surrounding_elements(1)
                                              CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,element_idx, &
                                                & INTERPOLATION_PARAMETERS(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
                                              en_idx=0
                                              XI_COORDINATES=0.0_DP
                                              number_of_nodes_xic(1)=DOMAIN%topology%elements%elements(element_idx)% &
                                                & basis%number_of_nodes_xic(1)
                                              number_of_nodes_xic(2)=DOMAIN%topology%elements%elements(element_idx)% &
                                                & basis%number_of_nodes_xic(2)
                                              IF(NUMBER_OF_DIMENSIONS==3) THEN
                                                number_of_nodes_xic(3)=DOMAIN%topology%elements%elements(element_idx)%basis% &
                                                  & number_of_nodes_xic(3)
                                              ELSE
                                                number_of_nodes_xic(3)=1
                                              ENDIF
!\todo: change definitions as soon as adjacent elements / boundary elements calculation works for simplex
                                              IF(DOMAIN%topology%elements%maximum_number_of_element_parameters==4.OR. &
                                                & DOMAIN%topology%elements%maximum_number_of_element_parameters==9.OR. &
                                                & DOMAIN%topology%elements%maximum_number_of_element_parameters==16.OR. &
                                                & DOMAIN%topology%elements%maximum_number_of_element_parameters==8.OR. &
                                                & DOMAIN%topology%elements%maximum_number_of_element_parameters==27.OR. &
                                                & DOMAIN%topology%elements%maximum_number_of_element_parameters==64) THEN
                                                  DO K=1,number_of_nodes_xic(3)
                                                    DO J=1,number_of_nodes_xic(2)
                                                      DO I=1,number_of_nodes_xic(1)
                                                        en_idx=en_idx+1
                                                          IF(DOMAIN%topology%elements%elements(element_idx)% &
                                                            & element_nodes(en_idx)==node_idx) EXIT
                                                          XI_COORDINATES(1)=XI_COORDINATES(1)+(1.0_DP/(number_of_nodes_xic(1)-1))
                                                      ENDDO
                                                      IF(DOMAIN%topology%elements%elements(element_idx)% &
                                                        & element_nodes(en_idx)==node_idx) EXIT
                                                        XI_COORDINATES(1)=0.0_DP
                                                        XI_COORDINATES(2)=XI_COORDINATES(2)+(1.0_DP/(number_of_nodes_xic(2)-1))
                                                    ENDDO
                                                    IF(DOMAIN%topology%elements%elements(element_idx)% &
                                                      & element_nodes(en_idx)==node_idx) EXIT
                                                    XI_COORDINATES(1)=0.0_DP
                                                    XI_COORDINATES(2)=0.0_DP
                                                    IF(number_of_nodes_xic(3)/=1) THEN
                                                      XI_COORDINATES(3)=XI_COORDINATES(3)+(1.0_DP/(number_of_nodes_xic(3)-1))
                                                    ENDIF
                                                  ENDDO
                                                  CALL FIELD_INTERPOLATE_XI(NO_PART_DERIV,XI_COORDINATES, &
                                                    & INTERPOLATED_POINT(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
                                              ELSE
!\todo: Use boundary flag
                                                IF(DOMAIN%topology%elements%maximum_number_of_element_parameters==3) THEN
                                                  T_COORDINATES(1,1:2)=[0.0_DP,1.0_DP]
                                                  T_COORDINATES(2,1:2)=[1.0_DP,0.0_DP]
                                                  T_COORDINATES(3,1:2)=[1.0_DP,1.0_DP]
                                                ELSE IF(DOMAIN%topology%elements%maximum_number_of_element_parameters==6) THEN
                                                  T_COORDINATES(1,1:2)=[0.0_DP,1.0_DP]
                                                  T_COORDINATES(2,1:2)=[1.0_DP,0.0_DP]
                                                  T_COORDINATES(3,1:2)=[1.0_DP,1.0_DP]
                                                  T_COORDINATES(4,1:2)=[0.5_DP,0.5_DP]
                                                  T_COORDINATES(5,1:2)=[1.0_DP,0.5_DP]
                                                  T_COORDINATES(6,1:2)=[0.5_DP,1.0_DP]
                                                ELSE IF(DOMAIN%topology%elements%maximum_number_of_element_parameters==10.AND. &
                                                  & NUMBER_OF_DIMENSIONS==2) THEN
                                                  T_COORDINATES(1,1:2)=[0.0_DP,1.0_DP]
                                                  T_COORDINATES(2,1:2)=[1.0_DP,0.0_DP]
                                                  T_COORDINATES(3,1:2)=[1.0_DP,1.0_DP]
                                                  T_COORDINATES(4,1:2)=[1.0_DP/3.0_DP,2.0_DP/3.0_DP]
                                                  T_COORDINATES(5,1:2)=[2.0_DP/3.0_DP,1.0_DP/3.0_DP]
                                                  T_COORDINATES(6,1:2)=[1.0_DP,1.0_DP/3.0_DP]
                                                  T_COORDINATES(7,1:2)=[1.0_DP,2.0_DP/3.0_DP]
                                                  T_COORDINATES(8,1:2)=[2.0_DP/3.0_DP,1.0_DP]
                                                  T_COORDINATES(9,1:2)=[1.0_DP/3.0_DP,1.0_DP]
                                                  T_COORDINATES(10,1:2)=[2.0_DP/3.0_DP,2.0_DP/3.0_DP]
                                                ELSE IF(DOMAIN%topology%elements%maximum_number_of_element_parameters==4) THEN
                                                  T_COORDINATES(1,1:3)=[0.0_DP,1.0_DP,1.0_DP]
                                                  T_COORDINATES(2,1:3)=[1.0_DP,0.0_DP,1.0_DP]
                                                  T_COORDINATES(3,1:3)=[1.0_DP,1.0_DP,0.0_DP]
                                                  T_COORDINATES(4,1:3)=[1.0_DP,1.0_DP,1.0_DP]
                                                ELSE IF(DOMAIN%topology%elements%maximum_number_of_element_parameters==10.AND. &
                                                  & NUMBER_OF_DIMENSIONS==3) THEN
                                                  T_COORDINATES(1,1:3)=[0.0_DP,1.0_DP,1.0_DP]
                                                  T_COORDINATES(2,1:3)=[1.0_DP,0.0_DP,1.0_DP]
                                                  T_COORDINATES(3,1:3)=[1.0_DP,1.0_DP,0.0_DP]
                                                  T_COORDINATES(4,1:3)=[1.0_DP,1.0_DP,1.0_DP]
                                                  T_COORDINATES(5,1:3)=[0.5_DP,0.5_DP,1.0_DP]
                                                  T_COORDINATES(6,1:3)=[0.5_DP,1.0_DP,0.5_DP]
                                                  T_COORDINATES(7,1:3)=[0.5_DP,1.0_DP,1.0_DP]
                                                  T_COORDINATES(8,1:3)=[1.0_DP,0.5_DP,0.5_DP]
                                                  T_COORDINATES(9,1:3)=[1.0_DP,1.0_DP,0.5_DP]
                                                  T_COORDINATES(10,1:3)=[1.0_DP,0.5_DP,1.0_DP]
                                                ELSE IF(DOMAIN%topology%elements%maximum_number_of_element_parameters==20) THEN
                                                  T_COORDINATES(1,1:3)=[0.0_DP,1.0_DP,1.0_DP]
                                                  T_COORDINATES(2,1:3)=[1.0_DP,0.0_DP,1.0_DP]
                                                  T_COORDINATES(3,1:3)=[1.0_DP,1.0_DP,0.0_DP]
                                                  T_COORDINATES(4,1:3)=[1.0_DP,1.0_DP,1.0_DP]
                                                  T_COORDINATES(5,1:3)=[1.0_DP/3.0_DP,2.0_DP/3.0_DP,1.0_DP]
                                                  T_COORDINATES(6,1:3)=[2.0_DP/3.0_DP,1.0_DP/3.0_DP,1.0_DP]
                                                  T_COORDINATES(7,1:3)=[1.0_DP/3.0_DP,1.0_DP,2.0_DP/3.0_DP]
                                                  T_COORDINATES(8,1:3)=[2.0_DP/3.0_DP,1.0_DP,1.0_DP/3.0_DP]
                                                  T_COORDINATES(9,1:3)=[1.0_DP/3.0_DP,1.0_DP,1.0_DP]
                                                  T_COORDINATES(10,1:3)=[2.0_DP/3.0_DP,1.0_DP,1.0_DP]
                                                  T_COORDINATES(11,1:3)=[1.0_DP,1.0_DP/3.0_DP,2.0_DP/3.0_DP]
                                                  T_COORDINATES(12,1:3)=[1.0_DP,2.0_DP/3.0_DP,1.0_DP/3.0_DP]
                                                  T_COORDINATES(13,1:3)=[1.0_DP,1.0_DP,1.0_DP/3.0_DP]
                                                  T_COORDINATES(14,1:3)=[1.0_DP,1.0_DP,2.0_DP/3.0_DP]
                                                  T_COORDINATES(15,1:3)=[1.0_DP,1.0_DP/3.0_DP,1.0_DP]
                                                  T_COORDINATES(16,1:3)=[1.0_DP,2.0_DP/3.0_DP,1.0_DP]
                                                  T_COORDINATES(17,1:3)=[2.0_DP/3.0_DP,2.0_DP/3.0_DP,2.0_DP/3.0_DP]
                                                  T_COORDINATES(18,1:3)=[2.0_DP/3.0_DP,2.0_DP/3.0_DP,1.0_DP]
                                                  T_COORDINATES(19,1:3)=[2.0_DP/3.0_DP,1.0_DP,2.0_DP/3.0_DP]
                                                  T_COORDINATES(20,1:3)=[1.0_DP,2.0_DP/3.0_DP,2.0_DP/3.0_DP]
                                                ENDIF
                                                DO K=1,DOMAIN%topology%elements%maximum_number_of_element_parameters
                                                  IF(DOMAIN%topology%elements%elements(element_idx)%element_nodes(K)==node_idx) EXIT
                                                ENDDO
                                                IF(NUMBER_OF_DIMENSIONS==2) THEN
                                                  CALL FIELD_INTERPOLATE_XI(NO_PART_DERIV,T_COORDINATES(K,1:2), &
                                                    & INTERPOLATED_POINT(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
                                                ELSE IF(NUMBER_OF_DIMENSIONS==3) THEN
                                                  CALL FIELD_INTERPOLATE_XI(NO_PART_DERIV,T_COORDINATES(K,1:3), &
                                                    & INTERPOLATED_POINT(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
                                                ENDIF
                                              ENDIF
                                              X=0.0_DP
                                              DO dim_idx=1,NUMBER_OF_DIMENSIONS
                                                X(dim_idx)=INTERPOLATED_POINT(FIELD_U_VARIABLE_TYPE)%ptr%VALUES(dim_idx,1)
                                              ENDDO !dim_idx
                                              !Loop over the derivatives
                                              CALL BOUNDARY_CONDITIONS_VARIABLE_GET(SOLVER_equations%BOUNDARY_CONDITIONS, &
                                                & dependentField%VARIABLE_TYPE_MAP(FIELD_U_VARIABLE_TYPE)%ptr, &
                                                & BOUNDARY_CONDITIONS_VARIABLE,err,error,*999)
                                              IF(ASSOCIATED(BOUNDARY_CONDITIONS_VARIABLE)) THEN
                                                DO deriv_idx=1,DOMAIN_NODES%NODES(node_idx)%NUMBER_OF_DERIVATIVES
                                                  ANALYTIC_FUNCTION_TYPE=EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE
                                                  GLOBAL_DERIV_INDEX=DOMAIN_NODES%NODES(node_idx)%DERIVATIVES(deriv_idx)% &
                                                    & GLOBAL_DERIVATIVE_INDEX
                                                  materialsField=>EQUATIONS_SET%MATERIALS%MATERIALS_FIELD
                                                  !Define MU_PARAM, density=1
                                                  MU_PARAM=materialsField%variables(1)%parameter_sets%parameter_sets(1)%ptr% &
                                                    & parameters%cmiss%data_dp(1)
                                                  !Define RHO_PARAM, density=2
                                                  IF(ANALYTIC_FUNCTION_TYPE==EQUATIONS_SET_STOKES_EQUATION_TWO_DIM_4.OR. &
                                                    & ANALYTIC_FUNCTION_TYPE==EQUATIONS_SET_STOKES_EQUATION_TWO_DIM_5.OR. &
                                                    & ANALYTIC_FUNCTION_TYPE==EQUATIONS_SET_STOKES_EQUATION_THREE_DIM_4.OR. &
                                                    & ANALYTIC_FUNCTION_TYPE==EQUATIONS_SET_STOKES_EQUATION_THREE_DIM_5) THEN
                                                    RHO_PARAM=materialsField%variables(1)%parameter_sets%parameter_sets(1)%ptr% &
                                                      & parameters%cmiss%data_dp(2)
                                                  ELSE
                                                    RHO_PARAM=0.0_DP
                                                  ENDIF
                                                  CALL STOKES_EQUATION_ANALYTIC_FUNCTIONS(VALUE,X,MU_PARAM,RHO_PARAM,CURRENT_TIME, &
                                                    & variable_type, &
                                                    & GLOBAL_DERIV_INDEX,ANALYTIC_FUNCTION_TYPE,NUMBER_OF_DIMENSIONS, &
                                                    & FIELD_VARIABLE%NUMBER_OF_COMPONENTS,component_idx,err,error,*999)
                                                  !Default to version 1 of each node derivative
                                                  local_ny=FIELD_VARIABLE%COMPONENTS(component_idx)%PARAM_TO_DOF_MAP% &
                                                    & NODE_PARAM2DOF_MAP%NODES(node_idx)%DERIVATIVES(deriv_idx)%VERSIONS(1)
                                                  CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(dependentField,variable_type, &
                                                    & FIELD_ANALYTIC_VALUES_SET_TYPE,local_ny,VALUE,err,error,*999)
                                                  BOUNDARY_CONDITION_CHECK_VARIABLE=BOUNDARY_CONDITIONS_VARIABLE% &
                                                    & CONDITION_TYPES(local_ny)
                                                  IF(BOUNDARY_CONDITION_CHECK_VARIABLE==BOUNDARY_CONDITION_FIXED) THEN
                                                   CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(dependentField, &
                                                     & variable_type,FIELD_VALUES_SET_TYPE,local_ny, &
                                                     & VALUE,err,error,*999)
                                                  ENDIF
                                                ENDDO !deriv_idx
                                              ELSE
                                                CALL FlagError("Boundary conditions U variable is not associated",err,error,*999)
                                              ENDIF
                                            ENDDO !node_idx
                                           ELSE
                                            CALL FlagError("Domain topology nodes is not associated.",err,error,*999)
                                          ENDIF
                                        ELSE
                                          CALL FlagError("Domain topology is not associated.",err,error,*999)
                                        ENDIF
                                      ELSE
                                        CALL FlagError("Domain is not associated.",err,error,*999)
                                      ENDIF
                                    ELSE
                                      CALL FlagError("Only node based interpolation is implemented.",err,error,*999)
                                    ENDIF
                                  ENDDO !component_idx
                                  CALL FIELD_PARAMETER_SET_UPDATE_START(dependentField,variable_type, &
                                   & FIELD_ANALYTIC_VALUES_SET_TYPE,err,error,*999)
                                  CALL FIELD_PARAMETER_SET_UPDATE_FINISH(dependentField,variable_type, &
                                   & FIELD_ANALYTIC_VALUES_SET_TYPE,err,error,*999)
                                  CALL FIELD_PARAMETER_SET_UPDATE_START(dependentField,variable_type, &
                                   & FIELD_VALUES_SET_TYPE,err,error,*999)
                                  CALL FIELD_PARAMETER_SET_UPDATE_FINISH(dependentField,variable_type, &
                                   & FIELD_VALUES_SET_TYPE,err,error,*999)
                                ELSE
                                  CALL FlagError("Field variable is not associated.",err,error,*999)
                                ENDIF
                               ENDDO !variable_idx
                               CALL FIELD_PARAMETER_SET_DATA_RESTORE(geometricField,FIELD_U_VARIABLE_TYPE,&
                                & FIELD_VALUES_SET_TYPE,GEOMETRIC_PARAMETERS,err,error,*999)
                            ELSE
                              CALL FlagError("Equations set geometric field is not associated.",err,error,*999)
                            ENDIF
                          ELSE
                            CALL FlagError("Equations set dependent field is not associated.",err,error,*999)
                          ENDIF
                        ELSE
                          CALL FlagError("Equations set analytic is not associated.",err,error,*999)
                        ENDIF
                      ELSE
                        CALL FlagError("Equations set is not associated.",err,error,*999)
                      ENDIF
                    ENDIF
                  ENDIF
                ELSE
                  CALL FlagError("Equations are not associated.",err,error,*999)
                END IF
              ELSE
                CALL FlagError("Solver equations are not associated.",err,error,*999)
              END IF
              CALL FIELD_PARAMETER_SET_UPDATE_START(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                & FIELD_VALUES_SET_TYPE,err,error,*999)
              CALL FIELD_PARAMETER_SET_UPDATE_FINISH(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                & FIELD_VALUES_SET_TYPE,err,error,*999)
            CASE(PROBLEM_PGM_STOKES_SUBTYPE)
             !Pre solve for the dynamic solver
             IF(SOLVER%SOLVE_TYPE==SOLVER_DYNAMIC_TYPE) THEN
               CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"Mesh movement change boundary conditions... ",err,error,*999)
                SOLVER_EQUATIONS=>SOLVER%SOLVER_EQUATIONS
                IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
                  SOLVER_MAPPING=>SOLVER_equations%SOLVER_MAPPING
                  EQUATIONS=>SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(1)%EQUATIONS
                  IF(ASSOCIATED(EQUATIONS)) THEN
                    EQUATIONS_SET=>equations%equationsSet
                    IF(ASSOCIATED(EQUATIONS_SET)) THEN
                      BOUNDARY_CONDITIONS=>SOLVER_equations%BOUNDARY_CONDITIONS
                      IF(ASSOCIATED(BOUNDARY_CONDITIONS)) THEN
                        CALL BOUNDARY_CONDITIONS_VARIABLE_GET(BOUNDARY_CONDITIONS,EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD% &
                          & VARIABLE_TYPE_MAP(FIELD_U_VARIABLE_TYPE)%ptr,BOUNDARY_CONDITIONS_VARIABLE,err,error,*999)
                        IF(ASSOCIATED(BOUNDARY_CONDITIONS_VARIABLE)) THEN
                          CALL FIELD_NUMBER_OF_COMPONENTS_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                            & NUMBER_OF_DIMENSIONS,err,error,*999)
                          NULLIFY(MESH_VELOCITY_VALUES)
                          CALL FIELD_PARAMETER_SET_DATA_GET(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                            & FIELD_MESH_VELOCITY_SET_TYPE,MESH_VELOCITY_VALUES,err,error,*999)
                          NULLIFY(BOUNDARY_VALUES)
                          CALL FIELD_PARAMETER_SET_DATA_GET(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                            & FIELD_BOUNDARY_SET_TYPE,BOUNDARY_VALUES,err,error,*999)
                          CALL FLUID_MECHANICS_IO_READ_BOUNDARY_CONDITIONS(SOLVER_LINEAR_TYPE,BOUNDARY_VALUES, &
                            & NUMBER_OF_DIMENSIONS,BOUNDARY_CONDITION_FIXED_INLET,CONTROL_LOOP%TIME_LOOP%INPUT_NUMBER, &
                            & CONTROL_LOOP%TIME_LOOP%ITERATION_NUMBER,CURRENT_TIME,1.0_DP,err,error,*999)
!                           DO equations_row_number=1,vectorEquations%vectorMapping%TOTAL_NUMBER_OF_ROWS
! xxxxxxxxxxxxxxxxxxxxxx
                          DO variable_idx=1,EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD%NUMBER_OF_VARIABLES
                            variable_type=EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD%VARIABLES(variable_idx)%variable_TYPE
                            FIELD_VARIABLE=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD%VARIABLE_TYPE_MAP(variable_type)%ptr
                            IF(ASSOCIATED(FIELD_VARIABLE)) THEN
                              DO component_idx=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
                                DOMAIN=>FIELD_VARIABLE%COMPONENTS(component_idx)%DOMAIN
                                IF(ASSOCIATED(DOMAIN)) THEN
                                  IF(ASSOCIATED(DOMAIN%TOPOLOGY)) THEN
                                    DOMAIN_NODES=>DOMAIN%TOPOLOGY%NODES
                                    IF(ASSOCIATED(DOMAIN_NODES)) THEN
                                      !Loop over the local nodes excluding the ghosts.
                                      DO node_idx=1,DOMAIN_NODES%NUMBER_OF_NODES
                                        DO deriv_idx=1,DOMAIN_NODES%NODES(node_idx)%NUMBER_OF_DERIVATIVES
                                          !Default to version 1 of each node derivative
                                          local_ny=FIELD_VARIABLE%COMPONENTS(component_idx)%PARAM_TO_DOF_MAP% &
                                            & NODE_PARAM2DOF_MAP%NODES(node_idx)%DERIVATIVES(deriv_idx)%VERSIONS(1)
! xxxxxxxxxxxxxxxxxxxxxxxxx
                                          DISPLACEMENT_VALUE=0.0_DP
                                          BOUNDARY_CONDITION_CHECK_VARIABLE=BOUNDARY_CONDITIONS_VARIABLE% &
                                            & CONDITION_TYPES(local_ny)
                                          IF(BOUNDARY_CONDITION_CHECK_VARIABLE==BOUNDARY_CONDITION_MOVED_WALL) THEN
                                            CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                                              & FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,local_ny, &
                                              & MESH_VELOCITY_VALUES(local_ny),err,error,*999)
                                          ELSE IF(BOUNDARY_CONDITION_CHECK_VARIABLE==BOUNDARY_CONDITION_FIXED_INLET) THEN
                                            CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                                              & FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,local_ny, &
                                              & BOUNDARY_VALUES(local_ny),err,error,*999)
                                          END IF
                                        ENDDO !deriv_idx
                                      ENDDO !node_idx
                                    ENDIF
                                  ENDIF
                                ENDIF
                              ENDDO !component_idx
                            ENDIF
                          ENDDO !variable_idx
                        ELSE
                          CALL FlagError("Boundary condition variable is not associated.",err,error,*999)
                        END IF
                      ELSE
                        CALL FlagError("Boundary conditions are not associated.",err,error,*999)
                      END IF
                    ELSE
                      CALL FlagError("Equations set is not associated.",err,error,*999)
                    END IF
                  ELSE
                    CALL FlagError("Equations are not associated.",err,error,*999)
                  END IF
                ELSE
                  CALL FlagError("Solver equations are not associated.",err,error,*999)
                END IF
                CALL FIELD_PARAMETER_SET_UPDATE_START(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_VALUES_SET_TYPE,err,error,*999)
                CALL FIELD_PARAMETER_SET_UPDATE_FINISH(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_VALUES_SET_TYPE,err,error,*999)
              END IF
            CASE(PROBLEM_ALE_STOKES_SUBTYPE)
              !Pre solve for the linear solver
              IF(SOLVER%SOLVE_TYPE==SOLVER_LINEAR_TYPE) THEN
                CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"Mesh movement change boundary conditions... ",err,error,*999)
                SOLVER_EQUATIONS=>SOLVER%SOLVER_EQUATIONS
                IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
                  SOLVER_MAPPING=>SOLVER_equations%SOLVER_MAPPING
                  EQUATIONS=>SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(1)%EQUATIONS
                  IF(ASSOCIATED(EQUATIONS)) THEN
                    EQUATIONS_SET=>equations%equationsSet
                    IF(ASSOCIATED(EQUATIONS_SET)) THEN
                      BOUNDARY_CONDITIONS=>SOLVER_equations%BOUNDARY_CONDITIONS
                      IF(ASSOCIATED(BOUNDARY_CONDITIONS)) THEN
                        CALL BOUNDARY_CONDITIONS_VARIABLE_GET(BOUNDARY_CONDITIONS,EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD% &
                          & VARIABLE_TYPE_MAP(FIELD_U_VARIABLE_TYPE)%ptr,BOUNDARY_CONDITIONS_VARIABLE,err,error,*999)
                        IF(ASSOCIATED(BOUNDARY_CONDITIONS_VARIABLE)) THEN
                          CALL FIELD_NUMBER_OF_COMPONENTS_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                            & NUMBER_OF_DIMENSIONS,err,error,*999)
                          NULLIFY(BOUNDARY_VALUES)
                          CALL FIELD_PARAMETER_SET_DATA_GET(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                            & FIELD_BOUNDARY_SET_TYPE,BOUNDARY_VALUES,err,error,*999)
                          CALL FLUID_MECHANICS_IO_READ_BOUNDARY_CONDITIONS(SOLVER_LINEAR_TYPE,BOUNDARY_VALUES, &
                            & NUMBER_OF_DIMENSIONS,BOUNDARY_CONDITION_MOVED_WALL,CONTROL_LOOP%TIME_LOOP%INPUT_NUMBER, &
                            & CONTROL_LOOP%TIME_LOOP%ITERATION_NUMBER,CURRENT_TIME,1.0_DP,err,error,*999)
                          DO variable_idx=1,EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD%NUMBER_OF_VARIABLES
                            variable_type=EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD%VARIABLES(variable_idx)%variable_TYPE
                            FIELD_VARIABLE=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD%VARIABLE_TYPE_MAP(variable_type)%ptr
                            IF(ASSOCIATED(FIELD_VARIABLE)) THEN
                              DO component_idx=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
                                DOMAIN=>FIELD_VARIABLE%COMPONENTS(component_idx)%DOMAIN
                                IF(ASSOCIATED(DOMAIN)) THEN
                                  IF(ASSOCIATED(DOMAIN%TOPOLOGY)) THEN
                                    DOMAIN_NODES=>DOMAIN%TOPOLOGY%NODES
                                    IF(ASSOCIATED(DOMAIN_NODES)) THEN
                                      !Loop over the local nodes excluding the ghosts.
                                      DO node_idx=1,DOMAIN_NODES%NUMBER_OF_NODES
                                        DO deriv_idx=1,DOMAIN_NODES%NODES(node_idx)%NUMBER_OF_DERIVATIVES
                                          !Default to version 1 of each node derivative
                                          local_ny=FIELD_VARIABLE%COMPONENTS(component_idx)%PARAM_TO_DOF_MAP% &
                                            & NODE_PARAM2DOF_MAP%NODES(node_idx)%DERIVATIVES(deriv_idx)%VERSIONS(1)
                                          BOUNDARY_CONDITION_CHECK_VARIABLE=BOUNDARY_CONDITIONS_VARIABLE% &
                                            & CONDITION_TYPES(local_ny)
                                          IF(BOUNDARY_CONDITION_CHECK_VARIABLE==BOUNDARY_CONDITION_MOVED_WALL) THEN
                                            CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                                              & FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,local_ny, &
                                              & BOUNDARY_VALUES(local_ny),err,error,*999)
                                          END IF
                                        END DO !deriv_idx
                                      ENDDO !node_idx
                                    ENDIF
                                  ENDIF
                                ENDIF
                              ENDDO !component_idx
                            ENDIF
                          ENDDO !variable_idx
                          CALL FIELD_PARAMETER_SET_DATA_RESTORE(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD, &
                            & FIELD_U_VARIABLE_TYPE,FIELD_BOUNDARY_SET_TYPE,BOUNDARY_VALUES,err,error,*999)
!\todo: This part should be read in out of a file eventually
                        ELSE
                          CALL FlagError("Boundary condition variable is not associated.",err,error,*999)
                        END IF
                      ELSE
                        CALL FlagError("Boundary conditions are not associated.",err,error,*999)
                      END IF
                    ELSE
                      CALL FlagError("Equations set is not associated.",err,error,*999)
                    END IF
                  ELSE
                    CALL FlagError("Equations are not associated.",err,error,*999)
                  END IF
                ELSE
                  CALL FlagError("Solver equations are not associated.",err,error,*999)
                END IF
                CALL FIELD_PARAMETER_SET_UPDATE_START(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_VALUES_SET_TYPE,err,error,*999)
                CALL FIELD_PARAMETER_SET_UPDATE_FINISH(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_VALUES_SET_TYPE,err,error,*999)
              !Pre solve for the dynamic solver
              ELSE IF(SOLVER%SOLVE_TYPE==SOLVER_DYNAMIC_TYPE) THEN
               CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"Mesh movement change boundary conditions... ",err,error,*999)
                SOLVER_EQUATIONS=>SOLVER%SOLVER_EQUATIONS
                IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
                  SOLVER_MAPPING=>SOLVER_equations%SOLVER_MAPPING
                  EQUATIONS=>SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(1)%EQUATIONS
                  IF(ASSOCIATED(EQUATIONS)) THEN
                    EQUATIONS_SET=>equations%equationsSet
                    IF(ASSOCIATED(EQUATIONS_SET)) THEN
                      BOUNDARY_CONDITIONS=>SOLVER_equations%BOUNDARY_CONDITIONS
                      IF(ASSOCIATED(BOUNDARY_CONDITIONS)) THEN
                        CALL BOUNDARY_CONDITIONS_VARIABLE_GET(BOUNDARY_CONDITIONS,EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD% &
                          & VARIABLE_TYPE_MAP(FIELD_U_VARIABLE_TYPE)%ptr,BOUNDARY_CONDITIONS_VARIABLE,err,error,*999)
                        IF(ASSOCIATED(BOUNDARY_CONDITIONS_VARIABLE)) THEN
                          CALL FIELD_NUMBER_OF_COMPONENTS_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                            & NUMBER_OF_DIMENSIONS,err,error,*999)
                          NULLIFY(MESH_VELOCITY_VALUES)
                          CALL FIELD_PARAMETER_SET_DATA_GET(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                            & FIELD_MESH_VELOCITY_SET_TYPE,MESH_VELOCITY_VALUES,err,error,*999)
                          NULLIFY(BOUNDARY_VALUES)
                          CALL FIELD_PARAMETER_SET_DATA_GET(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                            & FIELD_BOUNDARY_SET_TYPE,BOUNDARY_VALUES,err,error,*999)
                          CALL FLUID_MECHANICS_IO_READ_BOUNDARY_CONDITIONS(SOLVER_LINEAR_TYPE,BOUNDARY_VALUES, &
                            & NUMBER_OF_DIMENSIONS,BOUNDARY_CONDITION_FIXED_INLET,CONTROL_LOOP%TIME_LOOP%INPUT_NUMBER, &
                            & CONTROL_LOOP%TIME_LOOP%ITERATION_NUMBER,CURRENT_TIME,1.0_DP,err,error,*999)
                          DO variable_idx=1,EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD%NUMBER_OF_VARIABLES
                            variable_type=EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD%VARIABLES(variable_idx)%variable_TYPE
                            FIELD_VARIABLE=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD%VARIABLE_TYPE_MAP(variable_type)%ptr
                            IF(ASSOCIATED(FIELD_VARIABLE)) THEN
                              DO component_idx=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
                                DOMAIN=>FIELD_VARIABLE%COMPONENTS(component_idx)%DOMAIN
                                IF(ASSOCIATED(DOMAIN)) THEN
                                  IF(ASSOCIATED(DOMAIN%TOPOLOGY)) THEN
                                    DOMAIN_NODES=>DOMAIN%TOPOLOGY%NODES
                                    IF(ASSOCIATED(DOMAIN_NODES)) THEN
                                      !Loop over the local nodes excluding the ghosts.
                                      DO node_idx=1,DOMAIN_NODES%NUMBER_OF_NODES
                                        DO deriv_idx=1,DOMAIN_NODES%NODES(node_idx)%NUMBER_OF_DERIVATIVES
                                          !Default to version 1 of each node derivative
                                          local_ny=FIELD_VARIABLE%COMPONENTS(component_idx)%PARAM_TO_DOF_MAP% &
                                            & NODE_PARAM2DOF_MAP%NODES(node_idx)%DERIVATIVES(deriv_idx)%VERSIONS(1)
                                          DISPLACEMENT_VALUE=0.0_DP
                                          BOUNDARY_CONDITION_CHECK_VARIABLE=BOUNDARY_CONDITIONS_VARIABLE% &
                                            & CONDITION_TYPES(local_ny)
                                          IF(BOUNDARY_CONDITION_CHECK_VARIABLE==BOUNDARY_CONDITION_MOVED_WALL) THEN
                                            CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                                              & FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,local_ny, &
                                              & MESH_VELOCITY_VALUES(local_ny),err,error,*999)
                                          ELSE IF(BOUNDARY_CONDITION_CHECK_VARIABLE==BOUNDARY_CONDITION_FIXED_INLET) THEN
                                            CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                                              & FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,local_ny, &
                                              & BOUNDARY_VALUES(local_ny),err,error,*999)
                                          END IF
                                        END DO !deriv_idx
                                      ENDDO !node_idx
                                    ENDIF
                                  ENDIF
                                ENDIF
                              ENDDO !component_idx
                            ENDIF
                          ENDDO !variable_idx
                          CALL FIELD_PARAMETER_SET_DATA_RESTORE(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD, &
                            & FIELD_U_VARIABLE_TYPE,FIELD_MESH_VELOCITY_SET_TYPE,MESH_VELOCITY_VALUES,err,error,*999)
                          CALL FIELD_PARAMETER_SET_DATA_RESTORE(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD, &
                            & FIELD_U_VARIABLE_TYPE,FIELD_BOUNDARY_SET_TYPE,BOUNDARY_VALUES,err,error,*999)
                        ELSE
                          CALL FlagError("Boundary condition variable is not associated.",err,error,*999)
                        END IF
                      ELSE
                        CALL FlagError("Boundary conditions are not associated.",err,error,*999)
                      END IF
                    ELSE
                      CALL FlagError("Equations set is not associated.",err,error,*999)
                    END IF
                  ELSE
                    CALL FlagError("Equations are not associated.",err,error,*999)
                  END IF
                ELSE
                  CALL FlagError("Solver equations are not associated.",err,error,*999)
                END IF
                CALL FIELD_PARAMETER_SET_UPDATE_START(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_VALUES_SET_TYPE,err,error,*999)
                CALL FIELD_PARAMETER_SET_UPDATE_FINISH(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_VALUES_SET_TYPE,err,error,*999)
              END IF
              ! do nothing ???
            CASE DEFAULT
              localError="Problem subtype "//TRIM(NumberToVString(CONTROL_LOOP%PROBLEM%SPECIFICATION(3),"*",err,error))// &
                & " is not valid for a Stokes equation fluid type of a fluid mechanics problem class."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        ELSE
          CALL FlagError("Problem is not associated.",err,error,*999)
        ENDIF
      ELSE
        CALL FlagError("Solver is not associated.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Control loop is not associated.",err,error,*999)
    ENDIF
    EXITS("STOKES_PRE_SOLVE_UPDATE_BOUNDARY_CONDITIONS")
    RETURN
999 ERRORSEXITS("STOKES_PRE_SOLVE_UPDATE_BOUNDARY_CONDITIONS",err,error)
    RETURN 1
  END SUBROUTINE STOKES_PRE_SOLVE_UPDATE_BOUNDARY_CONDITIONS

  !
  !================================================================================================================================
  !
  !>Update mesh velocity and move mesh for ALE Stokes problem
  SUBROUTINE STOKES_PRE_SOLVE_ALE_UPDATE_MESH(CONTROL_LOOP,SOLVER,err,error,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP !<A pointer to the control loop to solve.
    TYPE(SOLVER_TYPE), POINTER :: SOLVER !<A pointer to the solvers
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(SOLVER_TYPE), POINTER :: SOLVER_ALE_STOKES, SOLVER_LAPLACE !<A pointer to the solvers
    TYPE(FIELD_TYPE), POINTER :: DEPENDENT_FIELD_LAPLACE, INDEPENDENT_FIELD_ALE_STOKES
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS_LAPLACE, SOLVER_EQUATIONS_ALE_STOKES  !<A pointer to the solver equations
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: SOLVER_MAPPING_LAPLACE, SOLVER_MAPPING_ALE_STOKES !<A pointer to the solver mapping
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET_LAPLACE, EQUATIONS_SET_ALE_STOKES !<A pointer to the equations set
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(VARYING_STRING) :: localError
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE
    TYPE(DOMAIN_TYPE), POINTER :: DOMAIN
    TYPE(DOMAIN_NODES_TYPE), POINTER :: DOMAIN_NODES

    REAL(DP) :: CURRENT_TIME,TIME_INCREMENT,ALPHA
    REAL(DP), POINTER :: MESH_DISPLACEMENT_VALUES(:)
    INTEGER(INTG) :: I,NUMBER_OF_DIMENSIONS_LAPLACE,NUMBER_OF_DIMENSIONS_ALE_STOKES,GEOMETRIC_MESH_COMPONENT
    INTEGER(INTG) :: INPUT_TYPE,INPUT_OPTION,component_idx,deriv_idx,local_ny,node_idx,variable_idx,variable_type

    ENTERS("STOKES_PRE_SOLVE_ALE_UPDATE_MESH",err,error,*999)

    IF(ASSOCIATED(CONTROL_LOOP)) THEN
      CALL CONTROL_LOOP_CURRENT_TIMES_GET(CONTROL_LOOP,CURRENT_TIME,TIME_INCREMENT,err,error,*999)
      NULLIFY(SOLVER_LAPLACE)
      NULLIFY(SOLVER_ALE_STOKES)
      IF(ASSOCIATED(SOLVER)) THEN
        IF(ASSOCIATED(CONTROL_LOOP%PROBLEM)) THEN
          IF(.NOT.ALLOCATED(control_loop%problem%specification)) THEN
            CALL FlagError("Problem specification is not allocated.",err,error,*999)
          ELSE IF(SIZE(control_loop%problem%specification,1)<3) THEN
            CALL FlagError("Problem specification must have three entries for a Stokes problem.",err,error,*999)
          END IF
          SELECT CASE(CONTROL_LOOP%PROBLEM%SPECIFICATION(3))
            CASE(PROBLEM_STATIC_STOKES_SUBTYPE,PROBLEM_LAPLACE_STOKES_SUBTYPE)
              ! do nothing ???
            CASE(PROBLEM_TRANSIENT_STOKES_SUBTYPE)
              ! do nothing ???
            CASE(PROBLEM_PGM_STOKES_SUBTYPE)
              !Update mesh within the dynamic solver
              IF(SOLVER%SOLVE_TYPE==SOLVER_DYNAMIC_TYPE) THEN
                !Get the independent field for the ALE Stokes problem
                CALL SOLVERS_SOLVER_GET(SOLVER%SOLVERS,1,SOLVER_ALE_STOKES,err,error,*999)
                SOLVER_EQUATIONS_ALE_STOKES=>SOLVER_ALE_STOKES%SOLVER_EQUATIONS
                IF(ASSOCIATED(SOLVER_EQUATIONS_ALE_STOKES)) THEN
                  SOLVER_MAPPING_ALE_STOKES=>SOLVER_EQUATIONS_ALE_STOKES%SOLVER_MAPPING
                  IF(ASSOCIATED(SOLVER_MAPPING_ALE_STOKES)) THEN
                    EQUATIONS_SET_ALE_STOKES=>SOLVER_MAPPING_ALE_STOKES%EQUATIONS_SETS(1)%ptr
                    IF(ASSOCIATED(EQUATIONS_SET_ALE_STOKES)) THEN
                      INDEPENDENT_FIELD_ALE_STOKES=>EQUATIONS_SET_ALE_STOKES%INDEPENDENT%INDEPENDENT_FIELD
                    ELSE
                      CALL FlagError("ALE Stokes equations set is not associated.",err,error,*999)
                    END IF
                    !Get the data
                    CALL FIELD_NUMBER_OF_COMPONENTS_GET(EQUATIONS_SET_ALE_STOKES%GEOMETRY%GEOMETRIC_FIELD, &
                      & FIELD_U_VARIABLE_TYPE,NUMBER_OF_DIMENSIONS_ALE_STOKES,err,error,*999)
!\todo: Introduce flags set by the user (42/1 only for testings purpose)
                    !Copy input to Stokes' independent field
                    INPUT_TYPE=42
                    INPUT_OPTION=1
                    NULLIFY(MESH_DISPLACEMENT_VALUES)
                    CALL FIELD_PARAMETER_SET_DATA_GET(EQUATIONS_SET_ALE_STOKES%INDEPENDENT%INDEPENDENT_FIELD, &
                      & FIELD_U_VARIABLE_TYPE,FIELD_MESH_DISPLACEMENT_SET_TYPE,MESH_DISPLACEMENT_VALUES,err,error,*999)
                    CALL FLUID_MECHANICS_IO_READ_DATA(SOLVER_LINEAR_TYPE,MESH_DISPLACEMENT_VALUES, &
                      & NUMBER_OF_DIMENSIONS_ALE_STOKES,INPUT_TYPE,INPUT_OPTION,CONTROL_LOOP%TIME_LOOP%ITERATION_NUMBER,1.0_DP, &
                      & err,error,*999)
                    CALL FIELD_PARAMETER_SET_UPDATE_START(EQUATIONS_SET_ALE_STOKES%INDEPENDENT%INDEPENDENT_FIELD, &
                      & FIELD_U_VARIABLE_TYPE,FIELD_MESH_DISPLACEMENT_SET_TYPE,err,error,*999)
                    CALL FIELD_PARAMETER_SET_UPDATE_FINISH(EQUATIONS_SET_ALE_STOKES%INDEPENDENT%INDEPENDENT_FIELD, &
                      & FIELD_U_VARIABLE_TYPE,FIELD_MESH_DISPLACEMENT_SET_TYPE,err,error,*999)
                  ELSE
                    CALL FlagError("ALE Stokes solver mapping is not associated.",err,error,*999)
                  END IF
                ELSE
                  CALL FlagError("ALE Stokes solver equations are not associated.",err,error,*999)
                END IF
                 !Use calculated values to update mesh
                CALL FIELD_COMPONENT_MESH_COMPONENT_GET(EQUATIONS_SET_ALE_STOKES%GEOMETRY%GEOMETRIC_FIELD, &
                  & FIELD_U_VARIABLE_TYPE,1,GEOMETRIC_MESH_COMPONENT,err,error,*999)
!                 CALL FIELD_PARAMETER_SET_DATA_GET(INDEPENDENT_FIELD_ALE_STOKES,FIELD_U_VARIABLE_TYPE, &
!                   & FIELD_MESH_DISPLACEMENT_SET_TYPE,MESH_DISPLACEMENT_VALUES,err,error,*999)
                EQUATIONS=>SOLVER_MAPPING_ALE_STOKES%EQUATIONS_SET_TO_SOLVER_MAP(1)%EQUATIONS
                IF(ASSOCIATED(EQUATIONS)) THEN
                  NULLIFY(vectorEquations)
                  CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
                  vectorMapping=>vectorEquations%vectorMapping
                  IF(ASSOCIATED(vectorMapping)) THEN
                    DO variable_idx=1,EQUATIONS_SET_ALE_STOKES%DEPENDENT%DEPENDENT_FIELD%NUMBER_OF_VARIABLES
                      variable_type=EQUATIONS_SET_ALE_STOKES%DEPENDENT%DEPENDENT_FIELD%VARIABLES(variable_idx)%variable_TYPE
                      FIELD_VARIABLE=>EQUATIONS_SET_ALE_STOKES%GEOMETRY%GEOMETRIC_FIELD%VARIABLE_TYPE_MAP(variable_type)%ptr
                      IF(ASSOCIATED(FIELD_VARIABLE)) THEN
                        DO component_idx=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
                          DOMAIN=>FIELD_VARIABLE%COMPONENTS(component_idx)%DOMAIN
                          IF(ASSOCIATED(DOMAIN)) THEN
                            IF(ASSOCIATED(DOMAIN%TOPOLOGY)) THEN
                              DOMAIN_NODES=>DOMAIN%TOPOLOGY%NODES
                              IF(ASSOCIATED(DOMAIN_NODES)) THEN
                                !Loop over the local nodes excluding the ghosts.
                                DO node_idx=1,DOMAIN_NODES%NUMBER_OF_NODES
                                  DO deriv_idx=1,DOMAIN_NODES%NODES(node_idx)%NUMBER_OF_DERIVATIVES
                                    !Default to version 1 of each node derivative
                                    local_ny=FIELD_VARIABLE%COMPONENTS(component_idx)%PARAM_TO_DOF_MAP% &
                                      & NODE_PARAM2DOF_MAP%NODES(node_idx)%DERIVATIVES(deriv_idx)%VERSIONS(1)
                                    CALL FIELD_PARAMETER_SET_ADD_LOCAL_DOF(EQUATIONS_SET_ALE_STOKES%GEOMETRY%GEOMETRIC_FIELD, &
                                      & FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,local_ny, &
                                      & MESH_DISPLACEMENT_VALUES(local_ny),err,error,*999)
                                  ENDDO !deriv_idx
                                ENDDO !node_idx
                              ENDIF
                            ENDIF
                          ENDIF
                        ENDDO !component_idx
                      ENDIF
                    ENDDO !variable_idx
                  ELSE
                    CALL FlagError("Equations mapping is not associated.",err,error,*999)
                  END IF
                ELSE
                  CALL FlagError("Equations are not associated.",err,error,*999)
                END IF
                CALL FIELD_PARAMETER_SET_UPDATE_START(EQUATIONS_SET_ALE_STOKES%GEOMETRY%GEOMETRIC_FIELD, &
                  & FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,err,error,*999)
                CALL FIELD_PARAMETER_SET_UPDATE_FINISH(EQUATIONS_SET_ALE_STOKES%GEOMETRY%GEOMETRIC_FIELD, &
                  & FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,err,error,*999)
                !Now use displacement values to calculate velocity values
                TIME_INCREMENT=CONTROL_LOOP%TIME_LOOP%TIME_INCREMENT
                ALPHA=1.0_DP/TIME_INCREMENT
                CALL FIELD_PARAMETER_SETS_COPY(INDEPENDENT_FIELD_ALE_STOKES,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_MESH_DISPLACEMENT_SET_TYPE,FIELD_MESH_VELOCITY_SET_TYPE,ALPHA,err,error,*999)
              ELSE
                CALL FlagError("Mesh motion calculation not successful for ALE problem.",err,error,*999)
              END IF
            CASE(PROBLEM_ALE_STOKES_SUBTYPE)
              !Update mesh within the dynamic solver
              IF(SOLVER%SOLVE_TYPE==SOLVER_DYNAMIC_TYPE) THEN
                IF(SOLVER%DYNAMIC_SOLVER%ALE) THEN
                  !Get the dependent field for the three component Laplace problem
                  CALL SOLVERS_SOLVER_GET(SOLVER%SOLVERS,1,SOLVER_LAPLACE,err,error,*999)
                  SOLVER_EQUATIONS_LAPLACE=>SOLVER_LAPLACE%SOLVER_EQUATIONS
                  IF(ASSOCIATED(SOLVER_EQUATIONS_LAPLACE)) THEN
                    SOLVER_MAPPING_LAPLACE=>SOLVER_EQUATIONS_LAPLACE%SOLVER_MAPPING
                    IF(ASSOCIATED(SOLVER_MAPPING_LAPLACE)) THEN
                      EQUATIONS_SET_LAPLACE=>SOLVER_MAPPING_LAPLACE%EQUATIONS_SETS(1)%ptr
                      IF(ASSOCIATED(EQUATIONS_SET_LAPLACE)) THEN
                        DEPENDENT_FIELD_LAPLACE=>EQUATIONS_SET_LAPLACE%DEPENDENT%DEPENDENT_FIELD
                      ELSE
                        CALL FlagError("Laplace equations set is not associated.",err,error,*999)
                      END IF
                      CALL FIELD_NUMBER_OF_COMPONENTS_GET(EQUATIONS_SET_LAPLACE%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                        & NUMBER_OF_DIMENSIONS_LAPLACE,err,error,*999)
                    ELSE
                      CALL FlagError("Laplace solver mapping is not associated.",err,error,*999)
                    END IF
                  ELSE
                    CALL FlagError("Laplace solver equations are not associated.",err,error,*999)
                  END IF
                  !Get the independent field for the ALE Stokes problem
                  CALL SOLVERS_SOLVER_GET(SOLVER%SOLVERS,2,SOLVER_ALE_STOKES,err,error,*999)
                  SOLVER_EQUATIONS_ALE_STOKES=>SOLVER_ALE_STOKES%SOLVER_EQUATIONS
                  IF(ASSOCIATED(SOLVER_EQUATIONS_ALE_STOKES)) THEN
                    SOLVER_MAPPING_ALE_STOKES=>SOLVER_EQUATIONS_ALE_STOKES%SOLVER_MAPPING
                    IF(ASSOCIATED(SOLVER_MAPPING_ALE_STOKES)) THEN
                      EQUATIONS_SET_ALE_STOKES=>SOLVER_MAPPING_ALE_STOKES%EQUATIONS_SETS(1)%ptr
                      IF(ASSOCIATED(EQUATIONS_SET_ALE_STOKES)) THEN
                        INDEPENDENT_FIELD_ALE_STOKES=>EQUATIONS_SET_ALE_STOKES%INDEPENDENT%INDEPENDENT_FIELD
                      ELSE
                        CALL FlagError("ALE Stokes equations set is not associated.",err,error,*999)
                      END IF
                      CALL FIELD_NUMBER_OF_COMPONENTS_GET(EQUATIONS_SET_ALE_STOKES%GEOMETRY%GEOMETRIC_FIELD, &
                        & FIELD_U_VARIABLE_TYPE,NUMBER_OF_DIMENSIONS_ALE_STOKES,err,error,*999)
                    ELSE
                      CALL FlagError("ALE Stokes solver mapping is not associated.",err,error,*999)
                    END IF
                  ELSE
                    CALL FlagError("ALE Stokes solver equations are not associated.",err,error,*999)
                  END IF
                  !Copy result from Laplace mesh movement to Stokes' independent field
                  IF(NUMBER_OF_DIMENSIONS_ALE_STOKES==NUMBER_OF_DIMENSIONS_LAPLACE) THEN
                    DO I=1,NUMBER_OF_DIMENSIONS_ALE_STOKES
                      CALL Field_ParametersToFieldParametersCopy(DEPENDENT_FIELD_LAPLACE, &
                        & FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,I,INDEPENDENT_FIELD_ALE_STOKES, &
                        & FIELD_U_VARIABLE_TYPE,FIELD_MESH_DISPLACEMENT_SET_TYPE,I,err,error,*999)
                    END DO
                  ELSE
                    CALL FlagError("Dimension of Laplace and ALE Stokes equations set is not consistent.",err,error,*999)
                  END IF
                  !Use calculated values to update mesh
                  CALL FIELD_COMPONENT_MESH_COMPONENT_GET(EQUATIONS_SET_ALE_STOKES%GEOMETRY%GEOMETRIC_FIELD, &
                    & FIELD_U_VARIABLE_TYPE,1,GEOMETRIC_MESH_COMPONENT,err,error,*999)
                  NULLIFY(MESH_DISPLACEMENT_VALUES)
                  CALL FIELD_PARAMETER_SET_DATA_GET(INDEPENDENT_FIELD_ALE_STOKES,FIELD_U_VARIABLE_TYPE, &
                    & FIELD_MESH_DISPLACEMENT_SET_TYPE,MESH_DISPLACEMENT_VALUES,err,error,*999)
                  EQUATIONS=>SOLVER_MAPPING_LAPLACE%EQUATIONS_SET_TO_SOLVER_MAP(1)%EQUATIONS
                  IF(ASSOCIATED(EQUATIONS)) THEN
                    NULLIFY(vectorEquations)
                    CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
                    vectorMapping=>vectorEquations%vectorMapping
                    IF(ASSOCIATED(vectorMapping)) THEN
                      DO variable_idx=1,EQUATIONS_SET_ALE_STOKES%DEPENDENT%DEPENDENT_FIELD%NUMBER_OF_VARIABLES
                        variable_type=EQUATIONS_SET_ALE_STOKES%DEPENDENT%DEPENDENT_FIELD%VARIABLES(variable_idx)%variable_TYPE
                        FIELD_VARIABLE=>EQUATIONS_SET_ALE_STOKES%GEOMETRY%GEOMETRIC_FIELD%VARIABLE_TYPE_MAP(variable_type)%ptr
                        IF(ASSOCIATED(FIELD_VARIABLE)) THEN
                          DO component_idx=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
                            DOMAIN=>FIELD_VARIABLE%COMPONENTS(component_idx)%DOMAIN
                            IF(ASSOCIATED(DOMAIN)) THEN
                              IF(ASSOCIATED(DOMAIN%TOPOLOGY)) THEN
                                DOMAIN_NODES=>DOMAIN%TOPOLOGY%NODES
                                IF(ASSOCIATED(DOMAIN_NODES)) THEN
                                  !Loop over the local nodes excluding the ghosts.
                                  DO node_idx=1,DOMAIN_NODES%NUMBER_OF_NODES
                                    DO deriv_idx=1,DOMAIN_NODES%NODES(node_idx)%NUMBER_OF_DERIVATIVES
                                      !Default to version 1 of each node derivative
                                      local_ny=FIELD_VARIABLE%COMPONENTS(component_idx)%PARAM_TO_DOF_MAP% &
                                        & NODE_PARAM2DOF_MAP%NODES(node_idx)%DERIVATIVES(deriv_idx)%VERSIONS(1)
                                      CALL FIELD_PARAMETER_SET_ADD_LOCAL_DOF(EQUATIONS_SET_ALE_STOKES%GEOMETRY%GEOMETRIC_FIELD, &
                                        & FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,local_ny, &
                                        & MESH_DISPLACEMENT_VALUES(local_ny),err,error,*999)
                                    ENDDO !deriv_idx
                                  ENDDO !node_idx
                                ENDIF
                              ENDIF
                            ENDIF
                          ENDDO !component_idx
                        ENDIF
                      ENDDO !variable_idx
                    ELSE
                      CALL FlagError("Equations mapping is not associated.",err,error,*999)
                    ENDIF
                    CALL FIELD_PARAMETER_SET_DATA_RESTORE(INDEPENDENT_FIELD_ALE_STOKES,FIELD_U_VARIABLE_TYPE, &
                      & FIELD_MESH_DISPLACEMENT_SET_TYPE,MESH_DISPLACEMENT_VALUES,err,error,*999)
                  ELSE
                    CALL FlagError("Equations are not associated.",err,error,*999)
                  END IF
                  CALL FIELD_PARAMETER_SET_UPDATE_START(EQUATIONS_SET_ALE_STOKES%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & FIELD_VALUES_SET_TYPE,err,error,*999)
                  CALL FIELD_PARAMETER_SET_UPDATE_FINISH(EQUATIONS_SET_ALE_STOKES%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & FIELD_VALUES_SET_TYPE,err,error,*999)
                  !Now use displacement values to calculate velocity values
                  TIME_INCREMENT=CONTROL_LOOP%TIME_LOOP%TIME_INCREMENT
                  ALPHA=1.0_DP/TIME_INCREMENT
                  CALL FIELD_PARAMETER_SETS_COPY(INDEPENDENT_FIELD_ALE_STOKES,FIELD_U_VARIABLE_TYPE, &
                    & FIELD_MESH_DISPLACEMENT_SET_TYPE,FIELD_MESH_VELOCITY_SET_TYPE,ALPHA,err,error,*999)
                ELSE
                  CALL FlagError("Mesh motion calculation not successful for ALE problem.",err,error,*999)
                END IF
              ELSE
                CALL FlagError("Mesh update is not defined for non-dynamic problems.",err,error,*999)
              END IF
            CASE DEFAULT
              localError="Problem subtype "//TRIM(NumberToVString(CONTROL_LOOP%PROBLEM%SPECIFICATION(3),"*",err,error))// &
                & " is not valid for a Stokes equation fluid type of a fluid mechanics problem class."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        ELSE
          CALL FlagError("Problem is not associated.",err,error,*999)
        ENDIF
      ELSE
        CALL FlagError("Solver is not associated.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Control loop is not associated.",err,error,*999)
    ENDIF
    EXITS("STOKES_PRE_SOLVE_ALE_UPDATE_MESH")
    RETURN
999 ERRORSEXITS("STOKES_PRE_SOLVE_ALE_UPDATE_MESH",err,error)
    RETURN 1
  END SUBROUTINE STOKES_PRE_SOLVE_ALE_UPDATE_MESH

  !
  !================================================================================================================================
  !
  !>Update mesh parameters for three component Laplace problem
  SUBROUTINE STOKES_PRE_SOLVE_ALE_UPDATE_PARAMETERS(CONTROL_LOOP,SOLVER,err,error,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP !<A pointer to the control loop to solve.
    TYPE(SOLVER_TYPE), POINTER :: SOLVER !<A pointer to the solver
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(FIELD_TYPE), POINTER :: independentField
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS  !<A pointer to the solver equations
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: SOLVER_MAPPING !<A pointer to the solver mapping
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set
    TYPE(EquationsType), POINTER :: EQUATIONS
    TYPE(VARYING_STRING) :: localError
    TYPE(DOMAIN_TYPE), POINTER :: DOMAIN
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE
    TYPE(DOMAIN_NODES_TYPE), POINTER :: DOMAIN_NODES

    REAL(DP) :: CURRENT_TIME,TIME_INCREMENT
    INTEGER(INTG) :: component_idx,node_idx,deriv_idx,local_ny,variable_idx,variable_type
    REAL(DP), POINTER :: MESH_STIFF_VALUES(:)


    ENTERS("STOKES_PRE_SOLVE_ALE_UPDATE_PARAMETERS",err,error,*999)

    IF(ASSOCIATED(CONTROL_LOOP)) THEN
      CALL CONTROL_LOOP_CURRENT_TIMES_GET(CONTROL_LOOP,CURRENT_TIME,TIME_INCREMENT,err,error,*999)
      IF(ASSOCIATED(SOLVER)) THEN
        IF(ASSOCIATED(CONTROL_LOOP%PROBLEM)) THEN
          IF(.NOT.ALLOCATED(control_loop%problem%specification)) THEN
            CALL FlagError("Problem specification is not allocated.",err,error,*999)
          ELSE IF(SIZE(control_loop%problem%specification,1)<3) THEN
            CALL FlagError("Problem specification must have three entries for a Stokes problem.",err,error,*999)
          END IF
          SELECT CASE(CONTROL_LOOP%PROBLEM%SPECIFICATION(3))
            CASE(PROBLEM_STATIC_STOKES_SUBTYPE,PROBLEM_LAPLACE_STOKES_SUBTYPE)
              ! do nothing ???
            CASE(PROBLEM_TRANSIENT_STOKES_SUBTYPE)
              ! do nothing ???
            CASE(PROBLEM_ALE_STOKES_SUBTYPE)
              IF(SOLVER%SOLVE_TYPE==SOLVER_LINEAR_TYPE) THEN
                !Get the independent field for the ALE Stokes problem
                SOLVER_EQUATIONS=>SOLVER%SOLVER_EQUATIONS
                IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
                  SOLVER_MAPPING=>SOLVER_equations%SOLVER_MAPPING
                  IF(ASSOCIATED(SOLVER_MAPPING)) THEN
                    EQUATIONS_SET=>SOLVER_MAPPING%EQUATIONS_SETS(1)%ptr
                    NULLIFY(MESH_STIFF_VALUES)
                    CALL FIELD_PARAMETER_SET_DATA_GET(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                      & FIELD_VALUES_SET_TYPE,MESH_STIFF_VALUES,err,error,*999)
                    IF(ASSOCIATED(EQUATIONS_SET)) THEN
                      EQUATIONS=>SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(1)%EQUATIONS
                      IF(ASSOCIATED(EQUATIONS)) THEN
                        independentField=>EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD
                        IF(ASSOCIATED(independentField)) THEN
                          DO variable_idx=1,EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD%NUMBER_OF_VARIABLES
                            variable_type=EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD%VARIABLES(variable_idx)%variable_TYPE
                            FIELD_VARIABLE=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD%VARIABLE_TYPE_MAP(variable_type)%ptr
                            IF(ASSOCIATED(FIELD_VARIABLE)) THEN
                              DO component_idx=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
                                DOMAIN=>FIELD_VARIABLE%COMPONENTS(component_idx)%DOMAIN
                                IF(ASSOCIATED(DOMAIN)) THEN
                                  IF(ASSOCIATED(DOMAIN%TOPOLOGY)) THEN
                                    DOMAIN_NODES=>DOMAIN%TOPOLOGY%NODES
                                    IF(ASSOCIATED(DOMAIN_NODES)) THEN
                                      !Loop over the local nodes excluding the ghosts.
                                      DO node_idx=1,DOMAIN_NODES%NUMBER_OF_NODES
                                        DO deriv_idx=1,DOMAIN_NODES%NODES(node_idx)%NUMBER_OF_DERIVATIVES
                                          !Default to version 1 of each node derivative
                                          local_ny=FIELD_VARIABLE%COMPONENTS(component_idx)%PARAM_TO_DOF_MAP% &
                                            & NODE_PARAM2DOF_MAP%NODES(node_idx)%DERIVATIVES(deriv_idx)%VERSIONS(1)
           !                             !Calculation of K values dependent on current mesh topology
                                          MESH_STIFF_VALUES(local_ny)=1.0_DP
                                          CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(EQUATIONS_SET%INDEPENDENT% &
                                            & INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,local_ny, &
                                            & MESH_STIFF_VALUES(local_ny),err,error,*999)
                                        ENDDO !deriv_idx
                                      ENDDO !node_idx
                                    ENDIF
                                  ENDIF
                                ENDIF
                              ENDDO !component_idx
                            ENDIF
                          ENDDO !variable_idx
                        ELSE
                          CALL FlagError("Independent field is not associated.",err,error,*999)
                        END IF
                      ELSE
                        CALL FlagError("Equations are not associated.",err,error,*999)
                      END IF
                    ELSE
                      CALL FlagError("Equations set is not associated.",err,error,*999)
                    ENDIF
                    CALL FIELD_PARAMETER_SET_DATA_RESTORE(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                      & FIELD_VALUES_SET_TYPE,MESH_STIFF_VALUES,err,error,*999)
                  ELSE
                    CALL FlagError("Solver mapping is not associated.",err,error,*999)
                  END IF
                ELSE
                  CALL FlagError("Solver equations are not associated.",err,error,*999)
                END IF
              ELSE IF(SOLVER%SOLVE_TYPE==SOLVER_DYNAMIC_TYPE) THEN
                CALL FlagError("Mesh motion calculation not successful for ALE problem.",err,error,*999)
              END IF
            CASE DEFAULT
              localError="Problem subtype "//TRIM(NumberToVString(CONTROL_LOOP%PROBLEM%SPECIFICATION(3),"*",err,error))// &
                & " is not valid for a Stokes equation fluid type of a fluid mechanics problem class."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        ELSE
          CALL FlagError("Problem is not associated.",err,error,*999)
        ENDIF
      ELSE
        CALL FlagError("Solver is not associated.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Control loop is not associated.",err,error,*999)
    ENDIF
    EXITS("STOKES_PRE_SOLVE_ALE_UPDATE_PARAMETERS")
    RETURN
999 ERRORSEXITS("STOKES_PRE_SOLVE_ALE_UPDATE_PARAMETERS",err,error)
    RETURN 1
  END SUBROUTINE STOKES_PRE_SOLVE_ALE_UPDATE_PARAMETERS

  !
  !================================================================================================================================
  !

  !>Output data post solve
  SUBROUTINE STOKES_POST_SOLVE_OUTPUT_DATA(CONTROL_LOOP,SOLVER,err,error,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP !<A pointer to the control loop to solve.
    TYPE(SOLVER_TYPE), POINTER :: SOLVER !<A pointer to the solver
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS  !<A pointer to the solver equations
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: SOLVER_MAPPING !<A pointer to the solver mapping
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set
    TYPE(FIELDS_TYPE), POINTER :: Fields
    TYPE(VARYING_STRING) :: localError,METHOD,FILENAME

    REAL(DP) :: CURRENT_TIME,TIME_INCREMENT
    INTEGER(INTG) :: EQUATIONS_SET_IDX,CURRENT_LOOP_ITERATION,OUTPUT_ITERATION_NUMBER,NUMBER_OF_DIMENSIONS
    LOGICAL :: EXPORT_FIELD
    CHARACTER(14) :: FILE,OUTPUT_FILE

    ENTERS("STOKES_POST_SOLVE_OUTPUT_DATA",err,error,*999)

    NULLIFY(SOLVER_EQUATIONS)
    NULLIFY(SOLVER_MAPPING)
    NULLIFY(EQUATIONS_SET)

    IF(ASSOCIATED(CONTROL_LOOP)) THEN
!       write(*,*)'CURRENT_TIME = ',CURRENT_TIME
!       write(*,*)'TIME_INCREMENT = ',TIME_INCREMENT
      IF(ASSOCIATED(SOLVER)) THEN
        IF(ASSOCIATED(CONTROL_LOOP%PROBLEM)) THEN
          IF(.NOT.ALLOCATED(control_loop%problem%specification)) THEN
            CALL FlagError("Problem specification is not allocated.",err,error,*999)
          ELSE IF(SIZE(control_loop%problem%specification,1)<3) THEN
            CALL FlagError("Problem specification must have three entries for a Stokes problem.",err,error,*999)
          END IF
          SELECT CASE(CONTROL_LOOP%PROBLEM%SPECIFICATION(3))
            CASE(PROBLEM_STATIC_STOKES_SUBTYPE,PROBLEM_LAPLACE_STOKES_SUBTYPE)
              SOLVER_EQUATIONS=>SOLVER%SOLVER_EQUATIONS
                IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
                  SOLVER_MAPPING=>SOLVER_equations%SOLVER_MAPPING
                  IF(ASSOCIATED(SOLVER_MAPPING)) THEN
                    !Make sure the equations sets are up to date
                    DO equations_set_idx=1,SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS
                      EQUATIONS_SET=>SOLVER_MAPPING%EQUATIONS_SETS(equations_set_idx)%ptr
                      METHOD="FORTRAN"
                      EXPORT_FIELD=.TRUE.
                      IF(EXPORT_FIELD) THEN
                        CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"...",err,error,*999)
                        CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"Now export fields... ",err,error,*999)
                        CALL FLUID_MECHANICS_IO_WRITE_CMGUI(EQUATIONS_SET%REGION,EQUATIONS_SET%GLOBAL_NUMBER,"STATICSOLUTION", &
                          & err,error,*999)
                        CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"STATICSOLUTION",err,error,*999)
                        CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"...",err,error,*999)
                      ENDIF
                    ENDDO
                  ENDIF
                ENDIF
            CASE(PROBLEM_TRANSIENT_STOKES_SUBTYPE,PROBLEM_ALE_STOKES_SUBTYPE,PROBLEM_PGM_STOKES_SUBTYPE)
              CALL CONTROL_LOOP_CURRENT_TIMES_GET(CONTROL_LOOP,CURRENT_TIME,TIME_INCREMENT,err,error,*999)
              SOLVER_EQUATIONS=>SOLVER%SOLVER_EQUATIONS
              IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
                SOLVER_MAPPING=>SOLVER_equations%SOLVER_MAPPING
                IF(ASSOCIATED(SOLVER_MAPPING)) THEN
                  !Make sure the equations sets are up to date
                  DO equations_set_idx=1,SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS
                    EQUATIONS_SET=>SOLVER_MAPPING%EQUATIONS_SETS(equations_set_idx)%ptr
                    CURRENT_LOOP_ITERATION=CONTROL_LOOP%TIME_LOOP%ITERATION_NUMBER
                    OUTPUT_ITERATION_NUMBER=CONTROL_LOOP%TIME_LOOP%OUTPUT_NUMBER
                    IF(OUTPUT_ITERATION_NUMBER/=0) THEN
                      IF(CONTROL_LOOP%TIME_LOOP%CURRENT_TIME<=CONTROL_LOOP%TIME_LOOP%STOP_TIME) THEN
                        IF(CURRENT_LOOP_ITERATION<10) THEN
                          WRITE(OUTPUT_FILE,'("TIME_STEP_000",I0)') CURRENT_LOOP_ITERATION
                        ELSE IF(CURRENT_LOOP_ITERATION<100) THEN
                          WRITE(OUTPUT_FILE,'("TIME_STEP_00",I0)') CURRENT_LOOP_ITERATION
                        ELSE IF(CURRENT_LOOP_ITERATION<1000) THEN
                          WRITE(OUTPUT_FILE,'("TIME_STEP_0",I0)') CURRENT_LOOP_ITERATION
                        ELSE IF(CURRENT_LOOP_ITERATION<10000) THEN
                          WRITE(OUTPUT_FILE,'("TIME_STEP_",I0)') CURRENT_LOOP_ITERATION
                        END IF

                        FILE=OUTPUT_FILE
                        FILENAME="./output/"//"MainTime_"//TRIM(NumberToVString(CURRENT_LOOP_ITERATION,"*",err,error))
                        METHOD="FORTRAN"
                        IF(MOD(CURRENT_LOOP_ITERATION,OUTPUT_ITERATION_NUMBER)==0)  THEN
                          IF(CONTROL_LOOP%outputtype >= CONTROL_LOOP_PROGRESS_OUTPUT) THEN
                            CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"...",err,error,*999)
                            CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"Now export fields... ",err,error,*999)
                          ENDIF
                          Fields=>EQUATIONS_SET%REGION%FIELDS
                          CALL FIELD_IO_NODES_EXPORT(Fields,FILENAME,METHOD,err,error,*999)
                          CALL FIELD_IO_ELEMENTS_EXPORT(Fields,FILENAME,METHOD,err,error,*999)
                          NULLIFY(Fields)
                          IF(CONTROL_LOOP%outputtype >= CONTROL_LOOP_PROGRESS_OUTPUT) THEN
                            CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,FILENAME,err,error,*999)
                            CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"...",err,error,*999)
                          ENDIF
                        END IF

                        IF(ASSOCIATED(EQUATIONS_SET%ANALYTIC)) THEN
                          IF(EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE==EQUATIONS_SET_NAVIER_STOKES_EQUATION_TWO_DIM_4.OR. &
                            & EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE==EQUATIONS_SET_NAVIER_STOKES_EQUATION_TWO_DIM_5.OR. &
                            & EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE==EQUATIONS_SET_NAVIER_STOKES_EQUATION_THREE_DIM_4.OR. &
                            & EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE==EQUATIONS_SET_NAVIER_STOKES_EQUATION_THREE_DIM_5.OR. &
                            & EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE==EQUATIONS_SET_NAVIER_STOKES_EQUATION_THREE_DIM_1) THEN
                            CALL AnalyticAnalysis_Output(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FILE,err,error,*999)
                          ENDIF
                        ENDIF
                      ENDIF
                    ENDIF
                  ENDDO
                ENDIF
              ENDIF
            CASE DEFAULT
              localError="Problem subtype "//TRIM(NumberToVString(CONTROL_LOOP%PROBLEM%SPECIFICATION(3),"*",err,error))// &
                & " is not valid for a Stokes equation fluid type of a fluid mechanics problem class."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        ELSE
          CALL FlagError("Problem is not associated.",err,error,*999)
        ENDIF
      ELSE
        CALL FlagError("Solver is not associated.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Control loop is not associated.",err,error,*999)
    ENDIF
    EXITS("STOKES_POST_SOLVE_OUTPUT_DATA")
    RETURN
999 ERRORSEXITS("STOKES_POST_SOLVE_OUTPUT_DATA",err,error)
    RETURN 1
  END SUBROUTINE STOKES_POST_SOLVE_OUTPUT_DATA

  !
  !================================================================================================================================
  !

  !>Calculates the analytic solution and sets the boundary conditions for an analytic problem.
  SUBROUTINE Stokes_BoundaryConditionsAnalyticCalculate(EQUATIONS_SET,BOUNDARY_CONDITIONS,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET
    TYPE(BOUNDARY_CONDITIONS_TYPE), POINTER :: BOUNDARY_CONDITIONS
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
!\todo: Reduce number of variables used
    INTEGER(INTG) :: component_idx,deriv_idx,dim_idx,local_ny,node_idx,NUMBER_OF_DIMENSIONS,variable_idx,variable_type,I,J,K
    INTEGER(INTG) :: number_of_nodes_xic(3),element_idx,en_idx,BOUND_COUNT,ANALYTIC_FUNCTION_TYPE,GLOBAL_DERIV_INDEX
    REAL(DP) :: VALUE,X(3),XI_COORDINATES(3)
!     REAL(DP) :: BOUNDARY_TOLERANCE, BOUNDARY_X(3,2),MU_PARAM,L
    REAL(DP) :: T_COORDINATES(20,3),CURRENT_TIME,MU_PARAM,RHO_PARAM
    REAL(DP), POINTER :: GEOMETRIC_PARAMETERS(:)
    TYPE(DOMAIN_TYPE), POINTER :: DOMAIN
    TYPE(DOMAIN_NODES_TYPE), POINTER :: DOMAIN_NODES
    TYPE(FIELD_TYPE), POINTER :: dependentField,geometricField,materialsField
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE,GEOMETRIC_VARIABLE
    TYPE(FIELD_INTERPOLATED_POINT_PTR_TYPE), POINTER :: INTERPOLATED_POINT(:)
    TYPE(FIELD_INTERPOLATION_PARAMETERS_PTR_TYPE), POINTER :: INTERPOLATION_PARAMETERS(:)
!     TYPE(VARYING_STRING) :: localError

! ! !     !Temp variables
! ! !     INTEGER(INTG) :: number_of_element_nodes,temp_local_ny,temp_node_number,velocity_DOF_check,temp_local_node_number

    ENTERS("Stokes_BoundaryConditionsAnalyticCalculate",err,error,*999)
!\todo: Introduce user call to set parameters
    BOUND_COUNT=0
! ! ! !     L=10.0_DP
    XI_COORDINATES(3)=0.0_DP
!     BOUNDARY_TOLERANCE=0.000000001_DP
! ! !     BOUNDARY_X=0.0_DP
! ! !     T_COORDINATES=0.0_DP
! ! !     number_of_element_nodes=0
! ! !     temp_local_node_number=0
! ! !     temp_local_ny=0
! ! !     temp_node_number=0
! ! !     velocity_DOF_check=0
    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(ASSOCIATED(EQUATIONS_SET%ANALYTIC)) THEN
        dependentField=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD
        IF(ASSOCIATED(dependentField)) THEN
          geometricField=>EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD
          IF(ASSOCIATED(geometricField)) THEN
            NULLIFY(INTERPOLATION_PARAMETERS)
            NULLIFY(INTERPOLATED_POINT)
            CALL FIELD_INTERPOLATION_PARAMETERS_INITIALISE(geometricField,INTERPOLATION_PARAMETERS,err,error,*999)
            CALL FIELD_INTERPOLATED_POINTS_INITIALISE(INTERPOLATION_PARAMETERS,INTERPOLATED_POINT,err,error,*999)
            CALL FIELD_NUMBER_OF_COMPONENTS_GET(geometricField,FIELD_U_VARIABLE_TYPE,NUMBER_OF_DIMENSIONS,err,error,*999)
! ! ! !\todo: Check adjacent element calculation / use boundary node flag instead / didn't work for simplex
! ! !             IF(NUMBER_OF_DIMENSIONS==2) THEN
! ! !               BOUNDARY_X(1,1)=0.0_DP
! ! !               BOUNDARY_X(1,2)=10.0_DP
! ! !               BOUNDARY_X(2,1)=0.0_DP
! ! !               BOUNDARY_X(2,2)=10.0_DP
! ! !             ELSE IF(NUMBER_OF_DIMENSIONS==3) THEN
! ! !               BOUNDARY_X(1,1)=-5.0_DP
! ! !               BOUNDARY_X(1,2)=5.0_DP
! ! !               BOUNDARY_X(2,1)=-5.0_DP
! ! !               BOUNDARY_X(2,2)=5.0_DP
! ! !               BOUNDARY_X(3,1)=-5.0_DP
! ! !               BOUNDARY_X(3,2)=5.0_DP
! ! !             ENDIF
            NULLIFY(GEOMETRIC_VARIABLE)
            CALL Field_VariableGet(geometricField,FIELD_U_VARIABLE_TYPE,GEOMETRIC_VARIABLE,err,error,*999)
            NULLIFY(GEOMETRIC_PARAMETERS)
            CALL FIELD_PARAMETER_SET_DATA_GET(geometricField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,GEOMETRIC_PARAMETERS, &
              & err,error,*999)
            IF(ASSOCIATED(BOUNDARY_CONDITIONS)) THEN
              DO variable_idx=1,dependentField%NUMBER_OF_VARIABLES
                variable_type=dependentField%VARIABLES(variable_idx)%variable_TYPE
                FIELD_VARIABLE=>dependentField%VARIABLE_TYPE_MAP(variable_type)%ptr
                IF(ASSOCIATED(FIELD_VARIABLE)) THEN
                  CALL FIELD_PARAMETER_SET_CREATE(dependentField,variable_type,FIELD_ANALYTIC_VALUES_SET_TYPE,err,error,*999)
                  DO component_idx=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
                    BOUND_COUNT=0
                    IF(FIELD_VARIABLE%COMPONENTS(component_idx)%INTERPOLATION_TYPE==FIELD_NODE_BASED_INTERPOLATION) THEN
                      DOMAIN=>FIELD_VARIABLE%COMPONENTS(component_idx)%DOMAIN
                      IF(ASSOCIATED(DOMAIN)) THEN
                        IF(ASSOCIATED(DOMAIN%TOPOLOGY)) THEN
                          DOMAIN_NODES=>DOMAIN%TOPOLOGY%NODES
                          IF(ASSOCIATED(DOMAIN_NODES)) THEN
                            !Loop over the local nodes excluding the ghosts.
                            DO node_idx=1,DOMAIN_NODES%NUMBER_OF_NODES
                              element_idx=DOMAIN%topology%nodes%nodes(node_idx)%surrounding_elements(1)
                              CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,element_idx, &
                                & INTERPOLATION_PARAMETERS(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
                              en_idx=0
                              XI_COORDINATES=0.0_DP
                              number_of_nodes_xic(1)=DOMAIN%topology%elements%elements(element_idx)%basis%number_of_nodes_xic(1)
                              number_of_nodes_xic(2)=DOMAIN%topology%elements%elements(element_idx)%basis%number_of_nodes_xic(2)
                              IF(NUMBER_OF_DIMENSIONS==3) THEN
                                number_of_nodes_xic(3)=DOMAIN%topology%elements%elements(element_idx)%basis%number_of_nodes_xic(3)
                              ELSE
                                number_of_nodes_xic(3)=1
                              ENDIF
  !\todo: Use boundary flag
                              IF(DOMAIN%topology%elements%maximum_number_of_element_parameters==4.AND.NUMBER_OF_DIMENSIONS==2 .OR. &
                                & DOMAIN%topology%elements%maximum_number_of_element_parameters==9.OR. &
                                & DOMAIN%topology%elements%maximum_number_of_element_parameters==16.OR. &
                                & DOMAIN%topology%elements%maximum_number_of_element_parameters==8.OR. &
                                & DOMAIN%topology%elements%maximum_number_of_element_parameters==27.OR. &
                                & DOMAIN%topology%elements%maximum_number_of_element_parameters==64) THEN
                                DO K=1,number_of_nodes_xic(3)
                                  DO J=1,number_of_nodes_xic(2)
                                    DO I=1,number_of_nodes_xic(1)
                                      en_idx=en_idx+1
                                      IF(DOMAIN%topology%elements%elements(element_idx)%element_nodes(en_idx)==node_idx) EXIT
                                      XI_COORDINATES(1)=XI_COORDINATES(1)+(1.0_DP/(number_of_nodes_xic(1)-1))
                                    ENDDO
                                      IF(DOMAIN%topology%elements%elements(element_idx)%element_nodes(en_idx)==node_idx) EXIT
                                      XI_COORDINATES(1)=0.0_DP
                                      XI_COORDINATES(2)=XI_COORDINATES(2)+(1.0_DP/(number_of_nodes_xic(2)-1))
                                  ENDDO
                                  IF(DOMAIN%topology%elements%elements(element_idx)%element_nodes(en_idx)==node_idx) EXIT
                                  XI_COORDINATES(1)=0.0_DP
                                  XI_COORDINATES(2)=0.0_DP
                                  IF(number_of_nodes_xic(3)/=1) THEN
                                    XI_COORDINATES(3)=XI_COORDINATES(3)+(1.0_DP/(number_of_nodes_xic(3)-1))
                                  ENDIF
                                ENDDO
                                CALL FIELD_INTERPOLATE_XI(NO_PART_DERIV,XI_COORDINATES, &
                                  & INTERPOLATED_POINT(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
                              ELSE
  !\todo: Use boundary flag
                                IF(DOMAIN%topology%elements%maximum_number_of_element_parameters==3) THEN
                                  T_COORDINATES(1,1:2)=[0.0_DP,1.0_DP]
                                  T_COORDINATES(2,1:2)=[1.0_DP,0.0_DP]
                                  T_COORDINATES(3,1:2)=[1.0_DP,1.0_DP]
                                ELSE IF(DOMAIN%topology%elements%maximum_number_of_element_parameters==6) THEN
                                  T_COORDINATES(1,1:2)=[0.0_DP,1.0_DP]
                                  T_COORDINATES(2,1:2)=[1.0_DP,0.0_DP]
                                  T_COORDINATES(3,1:2)=[1.0_DP,1.0_DP]
                                  T_COORDINATES(4,1:2)=[0.5_DP,0.5_DP]
                                  T_COORDINATES(5,1:2)=[1.0_DP,0.5_DP]
                                  T_COORDINATES(6,1:2)=[0.5_DP,1.0_DP]
                                ELSE IF(DOMAIN%topology%elements%maximum_number_of_element_parameters==10.AND. &
                                  & NUMBER_OF_DIMENSIONS==2) THEN
                                  T_COORDINATES(1,1:2)=[0.0_DP,1.0_DP]
                                  T_COORDINATES(2,1:2)=[1.0_DP,0.0_DP]
                                  T_COORDINATES(3,1:2)=[1.0_DP,1.0_DP]
                                  T_COORDINATES(4,1:2)=[1.0_DP/3.0_DP,2.0_DP/3.0_DP]
                                  T_COORDINATES(5,1:2)=[2.0_DP/3.0_DP,1.0_DP/3.0_DP]
                                  T_COORDINATES(6,1:2)=[1.0_DP,1.0_DP/3.0_DP]
                                  T_COORDINATES(7,1:2)=[1.0_DP,2.0_DP/3.0_DP]
                                  T_COORDINATES(8,1:2)=[2.0_DP/3.0_DP,1.0_DP]
                                  T_COORDINATES(9,1:2)=[1.0_DP/3.0_DP,1.0_DP]
                                  T_COORDINATES(10,1:2)=[2.0_DP/3.0_DP,2.0_DP/3.0_DP]
                                ELSE IF(DOMAIN%topology%elements%maximum_number_of_element_parameters==4) THEN
                                  T_COORDINATES(1,1:3)=[0.0_DP,1.0_DP,1.0_DP]
                                  T_COORDINATES(2,1:3)=[1.0_DP,0.0_DP,1.0_DP]
                                  T_COORDINATES(3,1:3)=[1.0_DP,1.0_DP,0.0_DP]
                                  T_COORDINATES(4,1:3)=[1.0_DP,1.0_DP,1.0_DP]
                                ELSE IF(DOMAIN%topology%elements%maximum_number_of_element_parameters==10.AND. &
                                  & NUMBER_OF_DIMENSIONS==3) THEN
                                  T_COORDINATES(1,1:3)=[0.0_DP,1.0_DP,1.0_DP]
                                  T_COORDINATES(2,1:3)=[1.0_DP,0.0_DP,1.0_DP]
                                  T_COORDINATES(3,1:3)=[1.0_DP,1.0_DP,0.0_DP]
                                  T_COORDINATES(4,1:3)=[1.0_DP,1.0_DP,1.0_DP]
                                  T_COORDINATES(5,1:3)=[0.5_DP,0.5_DP,1.0_DP]
                                  T_COORDINATES(6,1:3)=[0.5_DP,1.0_DP,0.5_DP]
                                  T_COORDINATES(7,1:3)=[0.5_DP,1.0_DP,1.0_DP]
                                  T_COORDINATES(8,1:3)=[1.0_DP,0.5_DP,0.5_DP]
                                  T_COORDINATES(9,1:3)=[1.0_DP,1.0_DP,0.5_DP]
                                  T_COORDINATES(10,1:3)=[1.0_DP,0.5_DP,1.0_DP]
                                ELSE IF(DOMAIN%topology%elements%maximum_number_of_element_parameters==20) THEN
                                  T_COORDINATES(1,1:3)=[0.0_DP,1.0_DP,1.0_DP]
                                  T_COORDINATES(2,1:3)=[1.0_DP,0.0_DP,1.0_DP]
                                  T_COORDINATES(3,1:3)=[1.0_DP,1.0_DP,0.0_DP]
                                  T_COORDINATES(4,1:3)=[1.0_DP,1.0_DP,1.0_DP]
                                  T_COORDINATES(5,1:3)=[1.0_DP/3.0_DP,2.0_DP/3.0_DP,1.0_DP]
                                  T_COORDINATES(6,1:3)=[2.0_DP/3.0_DP,1.0_DP/3.0_DP,1.0_DP]
                                  T_COORDINATES(7,1:3)=[1.0_DP/3.0_DP,1.0_DP,2.0_DP/3.0_DP]
                                  T_COORDINATES(8,1:3)=[2.0_DP/3.0_DP,1.0_DP,1.0_DP/3.0_DP]
                                  T_COORDINATES(9,1:3)=[1.0_DP/3.0_DP,1.0_DP,1.0_DP]
                                  T_COORDINATES(10,1:3)=[2.0_DP/3.0_DP,1.0_DP,1.0_DP]
                                  T_COORDINATES(11,1:3)=[1.0_DP,1.0_DP/3.0_DP,2.0_DP/3.0_DP]
                                  T_COORDINATES(12,1:3)=[1.0_DP,2.0_DP/3.0_DP,1.0_DP/3.0_DP]
                                  T_COORDINATES(13,1:3)=[1.0_DP,1.0_DP,1.0_DP/3.0_DP]
                                  T_COORDINATES(14,1:3)=[1.0_DP,1.0_DP,2.0_DP/3.0_DP]
                                  T_COORDINATES(15,1:3)=[1.0_DP,1.0_DP/3.0_DP,1.0_DP]
                                  T_COORDINATES(16,1:3)=[1.0_DP,2.0_DP/3.0_DP,1.0_DP]
                                  T_COORDINATES(17,1:3)=[2.0_DP/3.0_DP,2.0_DP/3.0_DP,2.0_DP/3.0_DP]
                                  T_COORDINATES(18,1:3)=[2.0_DP/3.0_DP,2.0_DP/3.0_DP,1.0_DP]
                                  T_COORDINATES(19,1:3)=[2.0_DP/3.0_DP,1.0_DP,2.0_DP/3.0_DP]
                                  T_COORDINATES(20,1:3)=[1.0_DP,2.0_DP/3.0_DP,2.0_DP/3.0_DP]
                                ENDIF
                                DO K=1,DOMAIN%topology%elements%maximum_number_of_element_parameters
                                  IF(DOMAIN%topology%elements%elements(element_idx)%element_nodes(K)==node_idx) EXIT
                                ENDDO
                                IF(NUMBER_OF_DIMENSIONS==2) THEN
                                  CALL FIELD_INTERPOLATE_XI(NO_PART_DERIV,T_COORDINATES(K,1:2), &
                                    & INTERPOLATED_POINT(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
                                ELSE IF(NUMBER_OF_DIMENSIONS==3) THEN
                                  CALL FIELD_INTERPOLATE_XI(NO_PART_DERIV,T_COORDINATES(K,1:3), &
                                    & INTERPOLATED_POINT(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
                                ENDIF
                              ENDIF
                              X=0.0_DP
                              DO dim_idx=1,NUMBER_OF_DIMENSIONS
                                X(dim_idx)=INTERPOLATED_POINT(FIELD_U_VARIABLE_TYPE)%ptr%VALUES(dim_idx,1)
                              ENDDO !dim_idx

                              !Loop over the derivatives
                              DO deriv_idx=1,DOMAIN_NODES%NODES(node_idx)%NUMBER_OF_DERIVATIVES
                                ANALYTIC_FUNCTION_TYPE=EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE
                                GLOBAL_DERIV_INDEX=DOMAIN_NODES%NODES(node_idx)%DERIVATIVES(deriv_idx)%GLOBAL_DERIVATIVE_INDEX
                                CURRENT_TIME=0.0_DP
                                materialsField=>EQUATIONS_SET%MATERIALS%MATERIALS_FIELD
                                !Define MU_PARAM, density=1
                                MU_PARAM=materialsField%variables(1)%parameter_sets%parameter_sets(1)%ptr% &
                                  & parameters%cmiss%data_dp(1)
                                !Define RHO_PARAM, density=2
                                IF(ANALYTIC_FUNCTION_TYPE==EQUATIONS_SET_STOKES_EQUATION_TWO_DIM_4.OR. &
                                  & ANALYTIC_FUNCTION_TYPE==EQUATIONS_SET_STOKES_EQUATION_TWO_DIM_5.OR. &
                                  & ANALYTIC_FUNCTION_TYPE==EQUATIONS_SET_STOKES_EQUATION_THREE_DIM_4.OR. &
                                  & ANALYTIC_FUNCTION_TYPE==EQUATIONS_SET_STOKES_EQUATION_THREE_DIM_5) THEN
                                  RHO_PARAM=materialsField%variables(1)%parameter_sets%parameter_sets(1)%ptr% &
                                    & parameters%cmiss%data_dp(2)
                                ELSE
                                  RHO_PARAM=0.0_DP
                                ENDIF
                                CALL STOKES_EQUATION_ANALYTIC_FUNCTIONS(VALUE,X,MU_PARAM,RHO_PARAM,CURRENT_TIME,variable_type, &
                                  & GLOBAL_DERIV_INDEX,ANALYTIC_FUNCTION_TYPE,NUMBER_OF_DIMENSIONS, &
                                  & FIELD_VARIABLE%NUMBER_OF_COMPONENTS,component_idx,err,error,*999)
                                !Default to version 1 of each node derivative
                                local_ny=FIELD_VARIABLE%COMPONENTS(component_idx)%PARAM_TO_DOF_MAP% &
                                  & NODE_PARAM2DOF_MAP%NODES(node_idx)%DERIVATIVES(deriv_idx)%VERSIONS(1)
                                CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(dependentField,variable_type, &
                                  & FIELD_ANALYTIC_VALUES_SET_TYPE,local_ny,VALUE,err,error,*999)
                                IF(variable_type==FIELD_U_VARIABLE_TYPE) THEN
  ! \todo: This part should work even for simplex elements as soon as adjacent element calculation has been fixed
                                  IF(DOMAIN_NODES%NODES(node_idx)%BOUNDARY_NODE) THEN
                                    !If we are a boundary node then set the analytic value on the boundary
                                    IF(component_idx<=NUMBER_OF_DIMENSIONS) THEN
                                      CALL BOUNDARY_CONDITIONS_SET_LOCAL_DOF(BOUNDARY_CONDITIONS,dependentField,variable_type, &
                                        & local_ny,BOUNDARY_CONDITION_FIXED,VALUE,err,error,*999)
                                      BOUND_COUNT=BOUND_COUNT+1
                                    ELSE
  ! \todo: This is just a workaround for linear pressure fields in simplex element components
                                      IF(DOMAIN%topology%elements%maximum_number_of_element_parameters==3) THEN
                                        IF(ANALYTIC_FUNCTION_TYPE==EQUATIONS_SET_STOKES_EQUATION_TWO_DIM_1.OR. &
                                          & ANALYTIC_FUNCTION_TYPE==EQUATIONS_SET_STOKES_EQUATION_TWO_DIM_2.OR. &
                                          & ANALYTIC_FUNCTION_TYPE==EQUATIONS_SET_STOKES_EQUATION_TWO_DIM_3.OR. &
                                          & ANALYTIC_FUNCTION_TYPE==EQUATIONS_SET_STOKES_EQUATION_TWO_DIM_4.OR. &
                                          & ANALYTIC_FUNCTION_TYPE==EQUATIONS_SET_STOKES_EQUATION_TWO_DIM_5) THEN
                                          IF(-0.001_DP<X(1).AND.X(1)<0.001_DP.AND.-0.001_DP<X(2).AND.X(2)<0.001_DP.OR. &
                                            &  10.0_DP-0.001_DP<X(1).AND.X(1)<10.0_DP+0.001_DP.AND.-0.001_DP<X(2).AND. &
                                            & X(2)<0.001_DP.OR. &
                                            &  10.0_DP-0.001_DP<X(1).AND.X(1)<10.0_DP+0.001_DP.AND.10.0_DP-0.001_DP<X(2).AND. &
                                            & X(2)<10.0_DP+0.001_DP.OR. &
                                            &  -0.001_DP<X(1).AND.X(1)<0.001_DP.AND.10.0_DP-0.001_DP<X(2).AND. &
                                            & X(2)<10.0_DP+0.001_DP) THEN
                                              CALL BOUNDARY_CONDITIONS_SET_LOCAL_DOF(BOUNDARY_CONDITIONS,dependentField, &
                                                & variable_type,local_ny,BOUNDARY_CONDITION_FIXED,VALUE,err,error,*999)
                                              BOUND_COUNT=BOUND_COUNT+1
                                          ENDIF
                                        ENDIF
                                      ELSE IF(DOMAIN%topology%elements%maximum_number_of_element_parameters==4.AND. &
                                        & NUMBER_OF_DIMENSIONS==3) THEN
                                        IF(ANALYTIC_FUNCTION_TYPE==EQUATIONS_SET_STOKES_EQUATION_THREE_DIM_1.OR. &
                                          & ANALYTIC_FUNCTION_TYPE==EQUATIONS_SET_STOKES_EQUATION_THREE_DIM_2.OR. &
                                          & ANALYTIC_FUNCTION_TYPE==EQUATIONS_SET_STOKES_EQUATION_THREE_DIM_3.OR. &
                                          & ANALYTIC_FUNCTION_TYPE==EQUATIONS_SET_STOKES_EQUATION_THREE_DIM_4.OR. &
                                          & ANALYTIC_FUNCTION_TYPE==EQUATIONS_SET_STOKES_EQUATION_THREE_DIM_5) THEN
                                          IF(-5.0_DP-0.001_DP<X(1).AND.X(1)<-5.0_DP+0.001_DP.AND.-5.0_DP-0.001_DP<X(2).AND. &
                                            & X(2)<-5.0_DP+0.001_DP.AND.-5.0_DP-0.001_DP<X(3).AND.X(3)<-5.0_DP+0.001_DP.OR. &
                                            & -5.0_DP-0.001_DP<X(1).AND.X(1)<-5.0_DP+0.001_DP.AND.5.0_DP-0.001_DP<X(2).AND. &
                                            & X(2)<5.0_DP+0.001_DP.AND.-5.0_DP-0.001_DP<X(3).AND.X(3)<-5.0_DP+0.001_DP.OR. &
                                            & 5.0_DP-0.001_DP<X(1).AND.X(1)<5.0_DP+0.001_DP.AND.5.0_DP-0.001_DP<X(2).AND. &
                                            & X(2)<5.0_DP+0.001_DP.AND.-5.0_DP-0.001_DP<X(3).AND.X(3)<-5.0_DP+0.001_DP.OR. &
                                            & 5.0_DP-0.001_DP<X(1).AND.X(1)<5.0_DP+0.001_DP.AND.-5.0_DP-0.001_DP<X(2).AND. &
                                            & X(2)<-5.0_DP+0.001_DP.AND.-5.0_DP-0.001_DP<X(3).AND.X(3)<-5.0_DP+0.001_DP.OR. &
                                            & -5.0_DP-0.001_DP<X(1).AND.X(1)<-5.0_DP+0.001_DP.AND.-5.0_DP-0.001_DP<X(2).AND. &
                                            & X(2)<-5.0_DP+0.001_DP.AND.5.0_DP-0.001_DP<X(3).AND.X(3)<5.0_DP+0.001_DP.OR. &
                                            & -5.0_DP-0.001_DP<X(1).AND.X(1)<-5.0_DP+0.001_DP.AND.5.0_DP-0.001_DP<X(2).AND. &
                                            & X(2)<5.0_DP+0.001_DP.AND.5.0_DP-0.001_DP<X(3).AND.X(3)<5.0_DP+0.001_DP.OR. &
                                            & 5.0_DP-0.001_DP<X(1).AND.X(1)<5.0_DP+0.001_DP.AND.5.0_DP-0.001_DP<X(2).AND. &
                                            & X(2)<5.0_DP+0.001_DP.AND.5.0_DP-0.001_DP<X(3).AND.X(3)<5.0_DP+0.001_DP.OR. &
                                            & 5.0_DP-0.001_DP<X(1).AND.X(1)<5.0_DP+0.001_DP.AND.-5.0_DP-0.001_DP<X(2).AND. &
                                            & X(2)<-5.0_DP+ 0.001_DP.AND.5.0_DP-0.001_DP<X(3).AND.X(3)<5.0_DP+0.001_DP) THEN
                                            CALL BOUNDARY_CONDITIONS_SET_LOCAL_DOF(BOUNDARY_CONDITIONS,dependentField, &
                                              & variable_type,local_ny,BOUNDARY_CONDITION_FIXED,VALUE,err,error,*999)
                                            BOUND_COUNT=BOUND_COUNT+1
                                          ENDIF
                                        ENDIF
  ! \todo: This is how it should be if adjacent elements would be working
                                      ELSE IF(BOUND_COUNT==0) THEN
                                        CALL BOUNDARY_CONDITIONS_SET_LOCAL_DOF(BOUNDARY_CONDITIONS,dependentField,variable_type, &
                                          & local_ny,BOUNDARY_CONDITION_FIXED,VALUE,err,error,*999)
                                        BOUND_COUNT=BOUND_COUNT+1
                                      ENDIF


                                    ENDIF
                                  ELSE
                                    IF(component_idx<=NUMBER_OF_DIMENSIONS) THEN
                                      CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(dependentField,variable_type, &
                                        & FIELD_VALUES_SET_TYPE,local_ny,VALUE,err,error,*999)
                                    ENDIF
                                  ENDIF
  ! \todo: Use boundary node flag
  ! ! !                                 !If we are a boundary node then set the analytic value on the boundary
  ! ! !                                 IF(NUMBER_OF_DIMENSIONS==2) THEN
  ! ! !                                   IF(X(1)<BOUNDARY_X(1,1)+BOUNDARY_TOLERANCE.AND.X(1)>BOUNDARY_X(1,1)-BOUNDARY_TOLERANCE.OR. &
  ! ! !                                     & X(1)<BOUNDARY_X(1,2)+BOUNDARY_TOLERANCE.AND.X(1)>BOUNDARY_X(1,2)-BOUNDARY_TOLERANCE.OR. &
  ! ! !                                     & X(2)<BOUNDARY_X(2,1)+BOUNDARY_TOLERANCE.AND.X(2)>BOUNDARY_X(2,1)-BOUNDARY_TOLERANCE.OR. &
  ! ! !                                     & X(2)<BOUNDARY_X(2,2)+BOUNDARY_TOLERANCE.AND.X(2)>BOUNDARY_X(2,2)-BOUNDARY_TOLERANCE) THEN
  ! ! !                                     IF(component_idx<=NUMBER_OF_DIMENSIONS) THEN
  ! ! !                                       CALL BOUNDARY_CONDITIONS_SET_LOCAL_DOF(BOUNDARY_CONDITIONS,variable_type,local_ny, &
  ! ! !                                         & BOUNDARY_CONDITION_FIXED,VALUE,err,error,*999)
  ! ! !                                     BOUND_COUNT=BOUND_COUNT+1
  ! ! !                                     !Apply boundary conditions check for pressure nodes
  ! ! !                                     ELSE IF(component_idx>NUMBER_OF_DIMENSIONS) THEN
  ! ! !                                       IF(DOMAIN%topology%elements%maximum_number_of_element_parameters==4) THEN
  ! ! !                                       IF(X(1)<BOUNDARY_X(1,1)+BOUNDARY_TOLERANCE.AND.X(1)>BOUNDARY_X(1,1)-BOUNDARY_TOLERANCE.AND. &
  ! ! !                                         & X(2)<BOUNDARY_X(2,1)+BOUNDARY_TOLERANCE.AND.X(2)>BOUNDARY_X(2,1)-BOUNDARY_TOLERANCE) &
  ! ! !                                         & THEN
  ! ! !                                            ! Commented out for testing purposes
  ! ! !                                           CALL BOUNDARY_CONDITIONS_SET_LOCAL_DOF(BOUNDARY_CONDITIONS,variable_type,local_ny, &
  ! ! !                                             & BOUNDARY_CONDITION_FIXED,VALUE,err,error,*999)
  ! ! !                                           BOUND_COUNT=BOUND_COUNT+1
  ! ! !                                       ENDIF
  ! ! !                                       ENDIF
  ! ! ! !\todo: Again, ...
  ! ! !                                       IF(DOMAIN%topology%elements%maximum_number_of_element_parameters==3.OR. &
  ! ! !                                         & DOMAIN%topology%elements%maximum_number_of_element_parameters==6.OR. &
  ! ! !                                         & DOMAIN%topology%elements%maximum_number_of_element_parameters==10) THEN
  ! ! !                                       IF(X(1)<BOUNDARY_X(1,1)+BOUNDARY_TOLERANCE.AND.X(1)>BOUNDARY_X(1,1)-BOUNDARY_TOLERANCE.AND. &
  ! ! !                                         & X(2)<BOUNDARY_X(2,1)+BOUNDARY_TOLERANCE.AND.X(2)>BOUNDARY_X(2,1)-BOUNDARY_TOLERANCE.OR. &
  ! ! !                                         & X(1)<BOUNDARY_X(1,1)+BOUNDARY_TOLERANCE.AND.X(1)>BOUNDARY_X(1,1)-BOUNDARY_TOLERANCE.AND.&
  ! ! !                                         & X(2)<BOUNDARY_X(2,2)+BOUNDARY_TOLERANCE.AND.X(2)>BOUNDARY_X(2,2)-BOUNDARY_TOLERANCE.OR. &
  ! ! !                                         & X(1)<BOUNDARY_X(1,2)+BOUNDARY_TOLERANCE.AND.X(1)>BOUNDARY_X(1,2)-BOUNDARY_TOLERANCE.AND.&
  ! ! !                                         & X(2)<BOUNDARY_X(2,1)+BOUNDARY_TOLERANCE.AND.X(2)>BOUNDARY_X(2,1)-BOUNDARY_TOLERANCE.OR. &
  ! ! !                                         & X(1)<BOUNDARY_X(1,2)+BOUNDARY_TOLERANCE.AND.X(1)>BOUNDARY_X(1,2)-BOUNDARY_TOLERANCE.AND.&
  ! ! !                                         & X(2)<BOUNDARY_X(2,2)+BOUNDARY_TOLERANCE.AND.X(2)>BOUNDARY_X(2,2)-BOUNDARY_TOLERANCE) &
  ! ! !                                         & THEN
  ! ! !                                           CALL BOUNDARY_CONDITIONS_SET_LOCAL_DOF(BOUNDARY_CONDITIONS,variable_type,local_ny, &
  ! ! !                                             & BOUNDARY_CONDITION_FIXED,VALUE,err,error,*999)
  ! ! !                                           BOUND_COUNT=BOUND_COUNT+1
  ! ! !                                       ENDIF
  ! ! !                                       ENDIF
  ! ! !                                     ENDIF
  ! ! !                                   ENDIF
  ! ! !                                     IF(component_idx<=NUMBER_OF_DIMENSIONS+1) THEN
  ! ! !                                       CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(dependentField,variable_type, &
  ! ! !                                         & FIELD_VALUES_SET_TYPE,local_ny,VALUE,err,error,*999)
  ! ! !                                     ENDIF
  ! ! !                                 ELSE IF(NUMBER_OF_DIMENSIONS==3) THEN
  ! ! !                                   IF(X(1)<BOUNDARY_X(1,1)+BOUNDARY_TOLERANCE.AND.X(1)>BOUNDARY_X(1,1)-BOUNDARY_TOLERANCE.OR. &
  ! ! !                                     & X(1)<BOUNDARY_X(1,2)+BOUNDARY_TOLERANCE.AND.X(1)>BOUNDARY_X(1,2)-BOUNDARY_TOLERANCE.OR. &
  ! ! !                                     & X(2)<BOUNDARY_X(2,1)+BOUNDARY_TOLERANCE.AND.X(2)>BOUNDARY_X(2,1)-BOUNDARY_TOLERANCE.OR. &
  ! ! !                                     & X(2)<BOUNDARY_X(2,2)+BOUNDARY_TOLERANCE.AND.X(2)>BOUNDARY_X(2,2)-BOUNDARY_TOLERANCE.OR. &
  ! ! !                                     & X(3)<BOUNDARY_X(3,1)+BOUNDARY_TOLERANCE.AND.X(3)>BOUNDARY_X(3,1)-BOUNDARY_TOLERANCE.OR. &
  ! ! !                                     & X(3)<BOUNDARY_X(3,2)+BOUNDARY_TOLERANCE.AND.X(3)>BOUNDARY_X(3,2)-BOUNDARY_TOLERANCE) THEN
  ! ! !                                     IF(component_idx<=NUMBER_OF_DIMENSIONS) THEN
  ! ! !                                       CALL BOUNDARY_CONDITIONS_SET_LOCAL_DOF(BOUNDARY_CONDITIONS,variable_type,local_ny, &
  ! ! !                                         & BOUNDARY_CONDITION_FIXED,VALUE,err,error,*999)
  ! ! !                                     BOUND_COUNT=BOUND_COUNT+1
  ! ! !                                     !Apply boundary conditions check for pressure nodes
  ! ! !                                     ELSE IF(component_idx>NUMBER_OF_DIMENSIONS) THEN
  ! ! !                                       IF(DOMAIN%topology%elements%maximum_number_of_element_parameters==4.OR. &
  ! ! !                                         & DOMAIN%topology%elements%maximum_number_of_element_parameters==10.OR. &
  ! ! !                                         & DOMAIN%topology%elements%maximum_number_of_element_parameters==20) THEN
  ! ! !                                       IF(X(1)<BOUNDARY_X(1,1)+BOUNDARY_TOLERANCE.AND.X(1)>BOUNDARY_X(1,1)-BOUNDARY_TOLERANCE.AND. &
  ! ! !                                        & X(2)<BOUNDARY_X(2,1)+BOUNDARY_TOLERANCE.AND.X(2)>BOUNDARY_X(2,1)-BOUNDARY_TOLERANCE.AND. &
  ! ! !                                        & X(3)<BOUNDARY_X(3,1)+BOUNDARY_TOLERANCE.AND.X(3)>BOUNDARY_X(3,1)-BOUNDARY_TOLERANCE.OR. &
  ! ! !                                        & X(1)<BOUNDARY_X(1,1)+BOUNDARY_TOLERANCE.AND.X(1)>BOUNDARY_X(1,1)-BOUNDARY_TOLERANCE.AND. &
  ! ! !                                        & X(2)<BOUNDARY_X(2,1)+BOUNDARY_TOLERANCE.AND.X(2)>BOUNDARY_X(2,1)-BOUNDARY_TOLERANCE.AND. &
  ! ! !                                        & X(3)<BOUNDARY_X(3,2)+BOUNDARY_TOLERANCE.AND.X(3)>BOUNDARY_X(3,2)-BOUNDARY_TOLERANCE.OR. &
  ! ! !                                        & X(1)<BOUNDARY_X(1,1)+BOUNDARY_TOLERANCE.AND.X(1)>BOUNDARY_X(1,1)-BOUNDARY_TOLERANCE.AND. &
  ! ! !                                        & X(2)<BOUNDARY_X(2,2)+BOUNDARY_TOLERANCE.AND.X(2)>BOUNDARY_X(2,2)-BOUNDARY_TOLERANCE.AND. &
  ! ! !                                        & X(3)<BOUNDARY_X(3,1)+BOUNDARY_TOLERANCE.AND.X(3)>BOUNDARY_X(3,1)-BOUNDARY_TOLERANCE.OR. &
  ! ! !                                        & X(1)<BOUNDARY_X(1,1)+BOUNDARY_TOLERANCE.AND.X(1)>BOUNDARY_X(1,1)-BOUNDARY_TOLERANCE.AND. &
  ! ! !                                        & X(2)<BOUNDARY_X(2,2)+BOUNDARY_TOLERANCE.AND.X(2)>BOUNDARY_X(2,2)-BOUNDARY_TOLERANCE.AND. &
  ! ! !                                        & X(3)<BOUNDARY_X(3,2)+BOUNDARY_TOLERANCE.AND.X(3)>BOUNDARY_X(3,2)-BOUNDARY_TOLERANCE.OR. &
  ! ! !                                        & X(1)<BOUNDARY_X(1,2)+BOUNDARY_TOLERANCE.AND.X(1)>BOUNDARY_X(1,2)-BOUNDARY_TOLERANCE.AND. &
  ! ! !                                        & X(2)<BOUNDARY_X(2,1)+BOUNDARY_TOLERANCE.AND.X(2)>BOUNDARY_X(2,1)-BOUNDARY_TOLERANCE.AND. &
  ! ! !                                        & X(3)<BOUNDARY_X(3,1)+BOUNDARY_TOLERANCE.AND.X(3)>BOUNDARY_X(3,1)-BOUNDARY_TOLERANCE.OR. &
  ! ! !                                        & X(1)<BOUNDARY_X(1,2)+BOUNDARY_TOLERANCE.AND.X(1)>BOUNDARY_X(1,2)-BOUNDARY_TOLERANCE.AND. &
  ! ! !                                        & X(2)<BOUNDARY_X(2,1)+BOUNDARY_TOLERANCE.AND.X(2)>BOUNDARY_X(2,1)-BOUNDARY_TOLERANCE.AND. &
  ! ! !                                        & X(3)<BOUNDARY_X(3,2)+BOUNDARY_TOLERANCE.AND.X(3)>BOUNDARY_X(3,2)-BOUNDARY_TOLERANCE.OR. &
  ! ! !                                        & X(1)<BOUNDARY_X(1,2)+BOUNDARY_TOLERANCE.AND.X(1)>BOUNDARY_X(1,2)-BOUNDARY_TOLERANCE.AND. &
  ! ! !                                        & X(2)<BOUNDARY_X(2,2)+BOUNDARY_TOLERANCE.AND.X(2)>BOUNDARY_X(2,2)-BOUNDARY_TOLERANCE.AND. &
  ! ! !                                        & X(3)<BOUNDARY_X(3,1)+BOUNDARY_TOLERANCE.AND.X(3)>BOUNDARY_X(3,1)-BOUNDARY_TOLERANCE.OR. &
  ! ! !                                        & X(1)<BOUNDARY_X(1,2)+BOUNDARY_TOLERANCE.AND.X(1)>BOUNDARY_X(1,2)-BOUNDARY_TOLERANCE.AND. &
  ! ! !                                        & X(2)<BOUNDARY_X(2,2)+BOUNDARY_TOLERANCE.AND.X(2)>BOUNDARY_X(2,2)-BOUNDARY_TOLERANCE.AND. &
  ! ! !                                        & X(3)<BOUNDARY_X(3,2)+BOUNDARY_TOLERANCE.AND.X(3)>BOUNDARY_X(3,2)-BOUNDARY_TOLERANCE) THEN
  ! ! !                                          CALL BOUNDARY_CONDITIONS_SET_LOCAL_DOF(BOUNDARY_CONDITIONS,variable_type,local_ny, &
  ! ! !                                            & BOUNDARY_CONDITION_FIXED,VALUE,err,error,*999)
  ! ! !                                          BOUND_COUNT=BOUND_COUNT+1
  ! ! !                                       ENDIF
  ! ! !                                       ENDIF
  ! ! !                                     ENDIF
  ! ! !                                   ELSE
  ! ! !                                     IF(component_idx<=NUMBER_OF_DIMENSIONS+1) THEN
  ! ! !                                       CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(dependentField,variable_type, &
  ! ! !                                         & FIELD_VALUES_SET_TYPE,local_ny,VALUE,err,error,*999)
  ! ! !                                     ENDIF
  ! ! !                                   ENDIF
  ! ! !                                 ENDIF
                                ENDIF
                              ENDDO !deriv_idx
                            ENDDO !node_idx
                          ELSE
                            CALL FlagError("Domain topology nodes is not associated.",err,error,*999)
                          ENDIF
                        ELSE
                          CALL FlagError("Domain topology is not associated.",err,error,*999)
                        ENDIF
                      ELSE
                        CALL FlagError("Domain is not associated.",err,error,*999)
                      ENDIF
                    ELSE
                      CALL FlagError("Only node based interpolation is implemented.",err,error,*999)
                    ENDIF
                  ENDDO !component_idx
                  CALL FIELD_PARAMETER_SET_UPDATE_START(dependentField,variable_type,FIELD_ANALYTIC_VALUES_SET_TYPE, &
                    & err,error,*999)
                  CALL FIELD_PARAMETER_SET_UPDATE_FINISH(dependentField,variable_type,FIELD_ANALYTIC_VALUES_SET_TYPE, &
                    & err,error,*999)
                  CALL FIELD_PARAMETER_SET_UPDATE_START(dependentField,variable_type,FIELD_VALUES_SET_TYPE, &
                    & err,error,*999)
                  CALL FIELD_PARAMETER_SET_UPDATE_FINISH(dependentField,variable_type,FIELD_VALUES_SET_TYPE, &
                    & err,error,*999)
                ELSE
                  CALL FlagError("Field variable is not associated.",err,error,*999)
                ENDIF
              ENDDO !variable_idx
              CALL FIELD_PARAMETER_SET_DATA_RESTORE(geometricField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                & GEOMETRIC_PARAMETERS,err,error,*999)
              CALL FIELD_INTERPOLATED_POINTS_FINALISE(INTERPOLATED_POINT,err,error,*999)
              CALL FIELD_INTERPOLATION_PARAMETERS_FINALISE(INTERPOLATION_PARAMETERS,err,error,*999)
            ELSE
              CALL FlagError("Boundary conditions is not associated.",err,error,*999)
            ENDIF
          ELSE
            CALL FlagError("Equations set geometric field is not associated.",err,error,*999)
          ENDIF
        ELSE
          CALL FlagError("Equations set dependent field is not associated.",err,error,*999)
        ENDIF
      ELSE
        CALL FlagError("Equations set analytic is not associated.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Equations set is not associated.",err,error,*999)
    ENDIF

    EXITS("Stokes_BoundaryConditionsAnalyticCalculate")
    RETURN
999 ERRORSEXITS("Stokes_BoundaryConditionsAnalyticCalculate",err,error)
    RETURN 1
  END SUBROUTINE Stokes_BoundaryConditionsAnalyticCalculate

  !
  !================================================================================================================================
  !
  !>Calculates the various analytic solutions given X and time, can be called from within analytic calculate or elsewhere if needed
  SUBROUTINE STOKES_EQUATION_ANALYTIC_FUNCTIONS(VALUE,X,MU_PARAM,RHO_PARAM,CURRENT_TIME,VARIABLE_TYPE, &
    & GLOBAL_DERIV_INDEX,ANALYTIC_FUNCTION_TYPE,NUMBER_OF_DIMENSIONS,NUMBER_OF_COMPONENTS,COMPONENT_IDX,err,error,*)

    !Argument variables
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    REAL(DP), INTENT(OUT) :: VALUE
    REAL(DP) :: MU_PARAM,RHO_PARAM
    REAL(DP), INTENT(IN) :: CURRENT_TIME
    REAL(DP), INTENT(IN), DIMENSION(3) :: X
    INTEGER(INTG), INTENT(IN) :: NUMBER_OF_DIMENSIONS,NUMBER_OF_COMPONENTS,COMPONENT_IDX
    !Local variables
    TYPE(VARYING_STRING) :: localError
    INTEGER(INTG) :: variable_type,GLOBAL_DERIV_INDEX,ANALYTIC_FUNCTION_TYPE
    !TYPE(DOMAIN_TYPE), POINTER :: DOMAIN
    !TYPE(DOMAIN_NODES_TYPE), POINTER :: DOMAIN_NODES
    REAL(DP) :: INTERNAL_TIME

    ENTERS("STOKES_EQUATION_ANALYTIC_FUNCTIONS",err,error,*999)

!\todo: Introduce user-defined or default values instead for density and viscosity
    INTERNAL_TIME=CURRENT_TIME
     SELECT CASE(ANALYTIC_FUNCTION_TYPE)
       CASE(EQUATIONS_SET_STOKES_EQUATION_TWO_DIM_1)
         IF(NUMBER_OF_DIMENSIONS==2.AND.NUMBER_OF_COMPONENTS==3) THEN
           !Polynomial function
           SELECT CASE(variable_type)
             CASE(FIELD_U_VARIABLE_TYPE)
               SELECT CASE(GLOBAL_DERIV_INDEX)
                 CASE(NO_GLOBAL_DERIV)
                   IF(component_idx==1) THEN
                     !calculate u
                     VALUE=X(2)**2/10.0_DP**2
                   ELSE IF(component_idx==2) THEN
                     !calculate v
                     VALUE=X(1)**2/10.0_DP**2
                   ELSE IF(component_idx==3) THEN
                     !calculate p
                     VALUE=2.0_DP*MU_PARAM/10.0_DP**2*X(1)
                   ELSE
                     CALL FlagError("Not implemented.",err,error,*999)
                   ENDIF
                 CASE(GLOBAL_DERIV_S1)
                   CALL FlagError("Not implemented.",err,error,*999)
                 CASE(GLOBAL_DERIV_S2)
                   CALL FlagError("Not implemented.",err,error,*999)
                 CASE(GLOBAL_DERIV_S1_S2)
                   CALL FlagError("Not implemented.",err,error,*999)
                 CASE DEFAULT
                   localError="The global derivative index of "//TRIM(NumberToVString( &
                     & GLOBAL_DERIV_INDEX,"*",err,error))// &
                     & " is invalid."
                   CALL FlagError(localError,err,error,*999)
               END SELECT
             CASE(FIELD_DELUDELN_VARIABLE_TYPE)
               SELECT CASE(GLOBAL_DERIV_INDEX)
                 CASE(NO_GLOBAL_DERIV)
                   VALUE= 0.0_DP
                 CASE(GLOBAL_DERIV_S1)
                   CALL FlagError("Not implemented.",err,error,*999)
                 CASE(GLOBAL_DERIV_S2)
                   CALL FlagError("Not implemented.",err,error,*999)
                 CASE(GLOBAL_DERIV_S1_S2)
                   CALL FlagError("Not implemented.",err,error,*999)
                 CASE DEFAULT
                   localError="The global derivative index of "//TRIM(NumberToVString( &
                     & GLOBAL_DERIV_INDEX,"*",err,error))// &
                     & " is invalid."
                   CALL FlagError(localError,err,error,*999)
               END SELECT
             CASE DEFAULT
               localError="The variable type of "//TRIM(NumberToVString(variable_type,"*",err,error))// &
                 & " is invalid."
               CALL FlagError(localError,err,error,*999)
           END SELECT
         ELSE
           localError="The number of components does not correspond to the number of dimensions."
           CALL FlagError(localError,err,error,*999)
         ENDIF
       CASE(EQUATIONS_SET_STOKES_EQUATION_TWO_DIM_2)
         IF(NUMBER_OF_DIMENSIONS==2.AND.NUMBER_OF_COMPONENTS==3) THEN
           !Exponential function
           SELECT CASE(variable_type)
             CASE(FIELD_U_VARIABLE_TYPE)
               SELECT CASE(GLOBAL_DERIV_INDEX)
                 CASE(NO_GLOBAL_DERIV)
                   IF(component_idx==1) THEN
                     !calculate u
                     VALUE= EXP((X(1)-X(2))/10.0_DP)
                   ELSE IF(component_idx==2) THEN
                     !calculate v
                     VALUE= EXP((X(1)-X(2))/10.0_DP)
                   ELSE IF(component_idx==3) THEN
                     !calculate p
                     VALUE= 2.0_DP*MU_PARAM/10.0_DP*EXP((X(1)-X(2))/10.0_DP)
                   ELSE
                     CALL FlagError("Not implemented.",err,error,*999)
                   ENDIF
                 CASE(GLOBAL_DERIV_S1)
                   CALL FlagError("Not implemented.",err,error,*999)
                 CASE(GLOBAL_DERIV_S2)
                   CALL FlagError("Not implemented.",err,error,*999)
                 CASE(GLOBAL_DERIV_S1_S2)
                   CALL FlagError("Not implemented.",err,error,*999)
                 CASE DEFAULT
                   localError="The global derivative index of "//TRIM(NumberToVString( &
                     & GLOBAL_DERIV_INDEX,"*",err,error))// &
                     & " is invalid."
                   CALL FlagError(localError,err,error,*999)
               END SELECT
             CASE(FIELD_DELUDELN_VARIABLE_TYPE)
               SELECT CASE(GLOBAL_DERIV_INDEX)
                 CASE(NO_GLOBAL_DERIV)
                   IF(component_idx==1) THEN
                     !calculate u
                     VALUE= 0.0_DP
                   ELSE IF(component_idx==2) THEN
                     !calculate v
                     VALUE= 0.0_DP
                   ELSE IF(component_idx==3) THEN
                     !calculate p
                     VALUE= 0.0_DP
                   ELSE
                     CALL FlagError("Not implemented.",err,error,*999)
                   ENDIF
                 CASE(GLOBAL_DERIV_S1)
                   CALL FlagError("Not implemented.",err,error,*999)
                 CASE(GLOBAL_DERIV_S2)
                   CALL FlagError("Not implemented.",err,error,*999)
                 CASE(GLOBAL_DERIV_S1_S2)
                   CALL FlagError("Not implemented.",err,error,*999)
                 CASE DEFAULT
                   localError="The global derivative index of "//TRIM(NumberToVString( &
                     & GLOBAL_DERIV_INDEX,"*",err,error))// &
                     & " is invalid."
                   CALL FlagError(localError,err,error,*999)
               END SELECT
             CASE DEFAULT
               localError="The variable type of "//TRIM(NumberToVString(variable_type,"*",err,error))// &
                 & " is invalid."
               CALL FlagError(localError,err,error,*999)
           END SELECT
         ELSE
           localError="The number of components does not correspond to the number of dimensions."
           CALL FlagError(localError,err,error,*999)
         ENDIF
       CASE(EQUATIONS_SET_STOKES_EQUATION_TWO_DIM_3)
         IF(NUMBER_OF_DIMENSIONS==2.AND.NUMBER_OF_COMPONENTS==3) THEN
           !Sine and cosine functions
           SELECT CASE(variable_type)
             CASE(FIELD_U_VARIABLE_TYPE)
               SELECT CASE(GLOBAL_DERIV_INDEX)
                 CASE(NO_GLOBAL_DERIV)
                   IF(component_idx==1) THEN
                     !calculate u
                     VALUE=SIN(2.0_DP*PI*X(1)/10.0_DP)*SIN(2.0_DP*PI*X(2)/10.0_DP)
                   ELSE IF(component_idx==2) THEN
                     !calculate v
                     VALUE=COS(2.0_DP*PI*X(1)/10.0_DP)*COS(2.0_DP*PI*X(2)/10.0_DP)
                   ELSE IF(component_idx==3) THEN
                     !calculate p
                     VALUE=4.0_DP*MU_PARAM*PI/10.0_DP*SIN(2.0_DP*PI*X(2)/10.0_DP)*COS(2.0_DP*PI*X(1)/10.0_DP)
                   ELSE
                     CALL FlagError("Not implemented.",err,error,*999)
                   ENDIF
                 CASE(GLOBAL_DERIV_S1)
                   CALL FlagError("Not implemented.",err,error,*999)
                 CASE(GLOBAL_DERIV_S2)
                   CALL FlagError("Not implemented.",err,error,*999)
                 CASE(GLOBAL_DERIV_S1_S2)
                   CALL FlagError("Not implemented.",err,error,*999)
                 CASE DEFAULT
                   localError="The global derivative index of "//TRIM(NumberToVString( &
                     & GLOBAL_DERIV_INDEX,"*",err,error))// &
                     & " is invalid."
                   CALL FlagError(localError,err,error,*999)
               END SELECT
             CASE(FIELD_DELUDELN_VARIABLE_TYPE)
               SELECT CASE(GLOBAL_DERIV_INDEX)
                 CASE(NO_GLOBAL_DERIV)
                   IF(component_idx==1) THEN
                     !calculate u
                     VALUE=0.0_DP
                   ELSE IF(component_idx==2) THEN
                     !calculate v
                     VALUE=16.0_DP*MU_PARAM*PI**2/10.0_DP**2*cos(2.0_DP*PI*X(2)/10.0_DP)*cos(2.0_DP*PI*X(1)/10.0_DP)
                   ELSE IF(component_idx==3) THEN
                     !calculate p
                     VALUE=0.0_DP
                   ELSE
                     CALL FlagError("Not implemented.",err,error,*999)
                   ENDIF
                 CASE(GLOBAL_DERIV_S1)
                   CALL FlagError("Not implemented.",err,error,*999)
                 CASE(GLOBAL_DERIV_S2)
                   CALL FlagError("Not implemented.",err,error,*999)
                 CASE(GLOBAL_DERIV_S1_S2)
                   CALL FlagError("Not implemented.",err,error,*999)
                 CASE DEFAULT
                   localError="The global derivative index of "//TRIM(NumberToVString( &
                     & GLOBAL_DERIV_INDEX,"*",err,error))// &
                     & " is invalid."
                   CALL FlagError(localError,err,error,*999)
               END SELECT
             CASE DEFAULT
               localError="The variable type of "//TRIM(NumberToVString(variable_type,"*",err,error))// &
                 & " is invalid."
               CALL FlagError(localError,err,error,*999)
           END SELECT
         ELSE
           localError="The number of components does not correspond to the number of dimensions."
           CALL FlagError(localError,err,error,*999)
         ENDIF
       CASE(EQUATIONS_SET_STOKES_EQUATION_TWO_DIM_4)
         IF(NUMBER_OF_DIMENSIONS==2.AND.NUMBER_OF_COMPONENTS==3) THEN
           !Reduced Taylor-Green solution for Stokes
           CALL FlagError("Not implemented.",err,error,*999)
         ENDIF
       CASE(EQUATIONS_SET_STOKES_EQUATION_TWO_DIM_5)
         IF(NUMBER_OF_DIMENSIONS==2.AND.NUMBER_OF_COMPONENTS==3) THEN
           !Stokes-Taylor-Green dynamic
           SELECT CASE(variable_type)
             CASE(FIELD_U_VARIABLE_TYPE)
               SELECT CASE(GLOBAL_DERIV_INDEX)
                 CASE(NO_GLOBAL_DERIV)
                   IF(component_idx==1) THEN
                     !calculate u
                     VALUE=X(2)*exp(-(2.0_DP*MU_PARAM/RHO_PARAM*CURRENT_TIME))
                   ELSE IF(component_idx==2) THEN
                     !calculate v
                     VALUE=X(1)*exp(-(2.0_DP*MU_PARAM/RHO_PARAM*CURRENT_TIME))
                   ELSE IF(component_idx==3) THEN
                     !calculate p
                     VALUE=2.0_DP*X(2)*MU_PARAM*exp(-(2.0_DP*MU_PARAM/RHO_PARAM*CURRENT_TIME))*X(1)
                   ELSE
                     CALL FlagError("Not implemented.",err,error,*999)
                   ENDIF
                 CASE(GLOBAL_DERIV_S1)
                   CALL FlagError("Not implemented.",err,error,*999)
                 CASE(GLOBAL_DERIV_S2)
                   CALL FlagError("Not implemented.",err,error,*999)
                 CASE(GLOBAL_DERIV_S1_S2)
                   CALL FlagError("Not implemented.",err,error,*999)
                 CASE DEFAULT
                   localError="The global derivative index of "//TRIM(NumberToVString( &
                     & GLOBAL_DERIV_INDEX,"*",err,error))// &
                     & " is invalid."
                   CALL FlagError(localError,err,error,*999)
               END SELECT
             CASE(FIELD_DELUDELN_VARIABLE_TYPE)
               SELECT CASE(GLOBAL_DERIV_INDEX)
                 CASE(NO_GLOBAL_DERIV)
                   IF(component_idx==1) THEN
                     !calculate u
                     VALUE=0.0_DP
                   ELSE IF(component_idx==2) THEN
                     !calculate v
                     VALUE=0.0_DP
                   ELSE IF(component_idx==3) THEN
                     !calculate p
                     VALUE=0.0_DP
                   ELSE
                     CALL FlagError("Not implemented.",err,error,*999)
                   ENDIF
                 CASE(GLOBAL_DERIV_S1)
                   CALL FlagError("Not implemented.",err,error,*999)
                 CASE(GLOBAL_DERIV_S2)
                   CALL FlagError("Not implemented.",err,error,*999)
                 CASE(GLOBAL_DERIV_S1_S2)
                   CALL FlagError("Not implemented.",err,error,*999)
                 CASE DEFAULT
                   localError="The global derivative index of "//TRIM(NumberToVString( &
                     & GLOBAL_DERIV_INDEX,"*",err,error))// &
                     & " is invalid."
                   CALL FlagError(localError,err,error,*999)
               END SELECT
             CASE DEFAULT
               localError="The variable type of "//TRIM(NumberToVString(variable_type,"*",err,error))// &
                 & " is invalid."
               CALL FlagError(localError,err,error,*999)
           END SELECT
         ELSE
           localError="The number of components does not correspond to the number of dimensions."
           CALL FlagError(localError,err,error,*999)
         ENDIF
       CASE(EQUATIONS_SET_STOKES_EQUATION_THREE_DIM_1)
         IF(NUMBER_OF_DIMENSIONS==3.AND.NUMBER_OF_COMPONENTS==4) THEN
!POLYNOM
           SELECT CASE(variable_type)
             CASE(FIELD_U_VARIABLE_TYPE)
               SELECT CASE(GLOBAL_DERIV_INDEX)
                 CASE(NO_GLOBAL_DERIV)
                   IF(component_idx==1) THEN
                     !calculate u
                     VALUE=X(2)**2/10.0_DP**2+X(3)**2/10.0_DP**2
                   ELSE IF(component_idx==2) THEN
                     !calculate v
                     VALUE=X(1)**2/10.0_DP**2+X(3)**2/10.0_DP**2
                   ELSE IF(component_idx==3) THEN
                     !calculate w
                     VALUE=X(1)**2/10.0_DP**2+X(2)**2/10.0_DP**2
                   ELSE IF(component_idx==4) THEN
                     !calculate p
                     VALUE=4.0_DP*MU_PARAM/10.0_DP**2*X(1)
                   ELSE
                     CALL FlagError("Not implemented.",err,error,*999)
                   ENDIF
                 CASE(GLOBAL_DERIV_S1)
                   CALL FlagError("Not implemented.",err,error,*999)
                 CASE(GLOBAL_DERIV_S2)
                   CALL FlagError("Not implemented.",err,error,*999)
                 CASE(GLOBAL_DERIV_S1_S2)
                   CALL FlagError("Not implemented.",err,error,*999)
                 CASE DEFAULT
                   localError="The global derivative index of "//TRIM(NumberToVString( &
                     & GLOBAL_DERIV_INDEX,"*",err,error))// &
                     & " is invalid."
                   CALL FlagError(localError,err,error,*999)
               END SELECT
             CASE(FIELD_DELUDELN_VARIABLE_TYPE)
               SELECT CASE(GLOBAL_DERIV_INDEX)
                 CASE(NO_GLOBAL_DERIV)
                   VALUE=0.0_DP
                 CASE(GLOBAL_DERIV_S1)
                   CALL FlagError("Not implemented.",err,error,*999)
                 CASE(GLOBAL_DERIV_S2)
                   CALL FlagError("Not implemented.",err,error,*999)
                 CASE(GLOBAL_DERIV_S1_S2)
                   CALL FlagError("Not implemented.",err,error,*999)
                 CASE DEFAULT
                   localError="The global derivative index of "//TRIM(NumberToVString( &
                     & GLOBAL_DERIV_INDEX,"*",err,error))// &
                     & " is invalid."
                   CALL FlagError(localError,err,error,*999)
                END SELECT
              CASE DEFAULT
                localError="The variable type of "//TRIM(NumberToVString(variable_type,"*",err,error))// &
                  & " is invalid."
                CALL FlagError(localError,err,error,*999)
           END SELECT
         ELSE
           localError="The number of components does not correspond to the number of dimensions."
           CALL FlagError(localError,err,error,*999)
         ENDIF
       CASE(EQUATIONS_SET_STOKES_EQUATION_THREE_DIM_2)
         IF(NUMBER_OF_DIMENSIONS==3.AND.NUMBER_OF_COMPONENTS==4) THEN
           !Exponential function
           SELECT CASE(variable_type)
             CASE(FIELD_U_VARIABLE_TYPE)
               SELECT CASE(GLOBAL_DERIV_INDEX)
                 CASE(NO_GLOBAL_DERIV)
                   IF(component_idx==1) THEN
                     !calculate u
                     VALUE=EXP((X(1)-X(2))/10.0_DP)+EXP((X(3)-X(1))/10.0_DP)
                   ELSE IF(component_idx==2) THEN
                     !calculate v
                     VALUE=EXP((X(1)-X(2))/10.0_DP)+EXP((X(2)-X(3))/10.0_DP)
                   ELSE IF(component_idx==3) THEN
                     !calculate w
                     VALUE=EXP((X(3)-X(1))/10.0_DP)+EXP((X(2)-X(3))/10.0_DP)
                   ELSE IF(component_idx==4) THEN
                     !calculate p
                     VALUE=2.0_DP*MU_PARAM/10.0_DP*(EXP((X(1)-X(2))/10.0_DP)-EXP((X(3)-X(1))/10.0_DP))
                   ELSE
                     CALL FlagError("Not implemented.",err,error,*999)
                   ENDIF
                 CASE(GLOBAL_DERIV_S1)
                   CALL FlagError("Not implemented.",err,error,*999)
                 CASE(GLOBAL_DERIV_S2)
                   CALL FlagError("Not implemented.",err,error,*999)
                 CASE(GLOBAL_DERIV_S1_S2)
                   CALL FlagError("Not implemented.",err,error,*999)
                 CASE DEFAULT
                   localError="The global derivative index of "//TRIM(NumberToVString( &
                     & GLOBAL_DERIV_INDEX,"*",err,error))// &
                     & " is invalid."
                   CALL FlagError(localError,err,error,*999)
               END SELECT
             CASE(FIELD_DELUDELN_VARIABLE_TYPE)
               SELECT CASE(GLOBAL_DERIV_INDEX)
                 CASE(NO_GLOBAL_DERIV)
                   IF(component_idx==1) THEN
                     !calculate u
                     VALUE=0.0_DP
                   ELSE IF(component_idx==2) THEN
                     !calculate v
                     VALUE=-2.0_DP*MU_PARAM*(2.0_DP*EXP(X(1)-X(2))+EXP(X(2)-X(3)))
                   ELSE IF(component_idx==3) THEN
                     !calculate w
                     VALUE=-2.0_DP*MU_PARAM*(2.0_DP*EXP(X(3)-X(1))+EXP(X(2)-X(3)))
                   ELSE IF(component_idx==4) THEN
                     !calculate p
                     VALUE=0.0_DP
                   ELSE
                     CALL FlagError("Not implemented.",err,error,*999)
                   ENDIF
                 CASE(GLOBAL_DERIV_S1)
                   CALL FlagError("Not implemented.",err,error,*999)
                 CASE(GLOBAL_DERIV_S2)
                   CALL FlagError("Not implemented.",err,error,*999)
                 CASE(GLOBAL_DERIV_S1_S2)
                   CALL FlagError("Not implemented.",err,error,*999)
                 CASE DEFAULT
                   localError="The global derivative index of "//TRIM(NumberToVString( &
                     & GLOBAL_DERIV_INDEX,"*",err,error))// &
                     & " is invalid."
                   CALL FlagError(localError,err,error,*999)
               END SELECT
             CASE DEFAULT
               localError="The variable type of "//TRIM(NumberToVString(variable_type,"*",err,error))// &
                 & " is invalid."
               CALL FlagError(localError,err,error,*999)
           END SELECT
         ELSE
           localError="The number of components does not correspond to the number of dimensions."
           CALL FlagError(localError,err,error,*999)
         ENDIF
       CASE(EQUATIONS_SET_STOKES_EQUATION_THREE_DIM_3)
         IF(NUMBER_OF_DIMENSIONS==3.AND.NUMBER_OF_COMPONENTS==4) THEN
           !Sine and cosine functions
           SELECT CASE(variable_type)
             CASE(FIELD_U_VARIABLE_TYPE)
               SELECT CASE(GLOBAL_DERIV_INDEX)
                 CASE(NO_GLOBAL_DERIV)
                   IF(component_idx==1) THEN
                     !calculate u
                     VALUE=sin(2.0_DP*PI*X(1)/10.0_DP)*sin(2.0_DP*PI*X(2)/10.0_DP)*sin(2.0_DP*PI*X(3)/10.0_DP)
                   ELSE IF(component_idx==2) THEN
                     !calculate v
                     VALUE=2.0_DP*cos(2.0_DP*PI*x(1)/10.0_DP)*sin(2.0_DP*PI*X(3)/10.0_DP)*cos(2.0_DP*PI*X(2)/10.0_DP)
                   ELSE IF(component_idx==3) THEN
                     !calculate w
                     VALUE=-cos(2.0_DP*PI*X(1)/10.0_DP)*sin(2.0_DP*PI*X(2)/10.0_DP)*cos(2.0_DP*PI*X(3)/10.0_DP)
                   ELSE IF(component_idx==4) THEN
                     !calculate p
                     VALUE=6.0_DP*MU_PARAM*PI/10.0_DP*sin(2.0_DP*PI*X(2)/10.0_DP)*sin(2.0_DP*PI*X(3)/10.0_DP)* &
                       & cos(2.0_DP*PI*X(1)/10.0_DP)
                   ELSE
                     CALL FlagError("Not implemented.",err,error,*999)
                   ENDIF
                 CASE(GLOBAL_DERIV_S1)
                   CALL FlagError("Not implemented.",err,error,*999)
                 CASE(GLOBAL_DERIV_S2)
                   CALL FlagError("Not implemented.",err,error,*999)
                 CASE(GLOBAL_DERIV_S1_S2)
                   CALL FlagError("Not implemented.",err,error,*999)
                 CASE DEFAULT
                   localError="The global derivative index of "//TRIM(NumberToVString( &
                     & GLOBAL_DERIV_INDEX,"*",err,error))// &
                     & " is invalid."
                   CALL FlagError(localError,err,error,*999)
                END SELECT
              CASE(FIELD_DELUDELN_VARIABLE_TYPE)
                SELECT CASE(GLOBAL_DERIV_INDEX)
                  CASE(NO_GLOBAL_DERIV)
                    IF(component_idx==1) THEN
                      !calculate u
                      VALUE=0.0_DP
                    ELSE IF(component_idx==2) THEN
                      !calculate v
                      VALUE=36*MU_PARAM*PI**2/10.0_DP**2*cos(2.0_DP*PI*X(2)/10.0_DP)*sin(2.0_DP*PI*X(3)/10.0_DP)* &
                        & cos(2.0_DP*PI*X(1)/10.0_DP)
                    ELSE IF(component_idx==3) THEN
                      !calculate w
                      VALUE=0.0_DP
                    ELSE IF(component_idx==4) THEN
                      !calculate p
                      VALUE=0.0_DP
                    ELSE
                      CALL FlagError("Not implemented.",err,error,*999)
                    ENDIF
                  CASE(GLOBAL_DERIV_S1)
                    CALL FlagError("Not implemented.",err,error,*999)
                  CASE(GLOBAL_DERIV_S2)
                    CALL FlagError("Not implemented.",err,error,*999)
                  CASE(GLOBAL_DERIV_S1_S2)
                    CALL FlagError("Not implemented.",err,error,*999)
                  CASE DEFAULT
                    localError="The global derivative index of "//TRIM(NumberToVString( &
                      & GLOBAL_DERIV_INDEX,"*",err,error))// &
                      & " is invalid."
                    CALL FlagError(localError,err,error,*999)
                END SELECT
              CASE DEFAULT
                localError="The variable type of "//TRIM(NumberToVString(variable_type,"*",err,error))// &
                  & " is invalid."
                CALL FlagError(localError,err,error,*999)
            END SELECT
          ELSE
            localError="The number of components does not correspond to the number of dimensions."
            CALL FlagError(localError,err,error,*999)
          ENDIF
       CASE(EQUATIONS_SET_STOKES_EQUATION_THREE_DIM_4)
         IF(NUMBER_OF_DIMENSIONS==3.AND.NUMBER_OF_COMPONENTS==4) THEN
           !Reduced Taylor-Green solution for Stokes
           CALL FlagError("Not implemented.",err,error,*999)
         ENDIF
       CASE(EQUATIONS_SET_STOKES_EQUATION_THREE_DIM_5)
         IF(NUMBER_OF_DIMENSIONS==3.AND.NUMBER_OF_COMPONENTS==4) THEN
           !Stokes-Taylor-Green dynamic
           SELECT CASE(variable_type)
             CASE(FIELD_U_VARIABLE_TYPE)
               SELECT CASE(GLOBAL_DERIV_INDEX)
                 CASE(NO_GLOBAL_DERIV)
                   IF(component_idx==1) THEN
                     !calculate u
                     VALUE=X(2)*exp(-(2.0_DP*MU_PARAM/RHO_PARAM*CURRENT_TIME))
                   ELSE IF(component_idx==2) THEN
                     !calculate v
                     VALUE=X(1)*exp(-(2.0_DP*MU_PARAM/RHO_PARAM*CURRENT_TIME))
                   ELSE IF(component_idx==3) THEN
                     !calculate v
                     VALUE=0.0_DP
                   ELSE IF(component_idx==4) THEN
                     !calculate p
                     VALUE=2.0_DP*X(2)*MU_PARAM*exp(-(2.0_DP*MU_PARAM/RHO_PARAM*CURRENT_TIME))*X(1)
                   ELSE
                     CALL FlagError("Not implemented.",err,error,*999)
                   ENDIF
                 CASE(GLOBAL_DERIV_S1)
                   CALL FlagError("Not implemented.",err,error,*999)
                 CASE(GLOBAL_DERIV_S2)
                   CALL FlagError("Not implemented.",err,error,*999)
                 CASE(GLOBAL_DERIV_S1_S2)
                   CALL FlagError("Not implemented.",err,error,*999)
                 CASE DEFAULT
                   localError="The global derivative index of "//TRIM(NumberToVString( &
                     & GLOBAL_DERIV_INDEX,"*",err,error))// &
                     & " is invalid."
                   CALL FlagError(localError,err,error,*999)
               END SELECT
             CASE(FIELD_DELUDELN_VARIABLE_TYPE)
               SELECT CASE(GLOBAL_DERIV_INDEX)
                 CASE(NO_GLOBAL_DERIV)
                   IF(component_idx==1) THEN
                     !calculate u
                     VALUE=0.0_DP
                   ELSE IF(component_idx==2) THEN
                     !calculate v
                     VALUE=0.0_DP
                   ELSE IF(component_idx==3) THEN
                     !calculate p
                     VALUE=0.0_DP
                   ELSE IF(component_idx==4) THEN
                     !calculate p
                     VALUE=0.0_DP
                   ELSE
                     CALL FlagError("Not implemented.",err,error,*999)
                   ENDIF
                 CASE(GLOBAL_DERIV_S1)
                   CALL FlagError("Not implemented.",err,error,*999)
                 CASE(GLOBAL_DERIV_S2)
                   CALL FlagError("Not implemented.",err,error,*999)
                 CASE(GLOBAL_DERIV_S1_S2)
                   CALL FlagError("Not implemented.",err,error,*999)
                 CASE DEFAULT
                   localError="The global derivative index of "//TRIM(NumberToVString( &
                     & GLOBAL_DERIV_INDEX,"*",err,error))// &
                     & " is invalid."
                   CALL FlagError(localError,err,error,*999)
               END SELECT
             CASE DEFAULT
               localError="The variable type of "//TRIM(NumberToVString(variable_type,"*",err,error))// &
                 & " is invalid."
               CALL FlagError(localError,err,error,*999)
           END SELECT
         ELSE
           localError="The number of components does not correspond to the number of dimensions."
           CALL FlagError(localError,err,error,*999)
         ENDIF
        CASE DEFAULT
          localError="The analytic function type of "// &
            & TRIM(NumberToVString(ANALYTIC_FUNCTION_TYPE,"*",err,error))// &
            & " is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT

    EXITS("STOKES_EQUATION_ANALYTIC_FUNCTIONS")
    RETURN
999 ERRORSEXITS("STOKES_EQUATION_ANALYTIC_FUNCTIONS",err,error)
    RETURN 1

  END SUBROUTINE STOKES_EQUATION_ANALYTIC_FUNCTIONS

  !
  !================================================================================================================================
  !

END MODULE STOKES_EQUATIONS_ROUTINES
