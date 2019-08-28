################################################################################
include /vol/bob/make/head.mk
################################################################################

BUILD := "$(shell perl -e '$$pwd=$$ENV{PWD}; $$time=localtime(); print "\\\"$$pwd $$time\\\"";')"

IPATHS 	+= -I$(MKLROOT)/include


################################################################################

LIB1       := Wstrandsim

LIB1_CF    += -DBUILD=${BUILD} -DWITH_MKL -DHAVE_OPENMP -DTIMING_ON
LIB1_CF    +=  -rdynamic


ifneq "$(filter LI%,$(ABI))" ""
LIB1_CF += -Werror
LIB1_CF    += -openmp -DHAVE_OPENMP
LIB1_CF    += -fast -axSSE4.2 # I don't know what the equivalent in g++ would be
#LIB1_CF += -wd21, # type qualifiers are meaningless in this declaration
LIB1_CF += -wd473 # returning reference to local temporary
LIB1_CF += -wd1125 # entity-kind "entity" is hidden by "entity" -- virtual function override intended?
LIB1_CF += -wd2259 # non-pointer conversion from "type" to "type" may lose significant bits
LIB1_CF += -wd2304 # non-explicit constructor with single argument may cause implicit type conversion
LIB1_CF += -wd2305 # declaration of 'explicit' constructor without a single argument is redundant
LIB1_CF += -Weffc++ # Turn on Effective C++ warnings
LIB1_CF += -wd2014 # Effective C++ Item 3 prefer new and delete to malloc and free
LIB1_CF += -wd2015 # Effective C++ Item 4 prefer C++ style comments
LIB1_CF += -wd2021 # Effective C++ Item 11 declare a copy constructor and assignment operator for xxxx
LIB1_CF += -wd2022 # Effective C++ Item 12 field entity-kind "entity" not initialized (preferable to assignment in constructors)
else
ifneq "$(filter L4%,$(filter-out L40,$(ABI)))" ""
LIB1_CF    += -fopenmp -DHAVE_OPENMP
endif
endif

# If we are not in release mode, we like to enable assertions even in opt
ifndef release
ifndef ndebug
LIB1_CF += -UNDEBUG -g
endif
endif

LIB1_HSD   := Wstrandsim

# Header files
LIB1_HDIR := .
LIB1_H += Collision/BoundingBox.hh
LIB1_H += Collision/BVH.hh
LIB1_H += Collision/BVHNode.hh
LIB1_H += Collision/CollisionBase.hh
LIB1_H += Collision/CollisionDetector.hh
LIB1_H += Collision/CollisionUtils.hh
LIB1_H += Collision/ContinuousTimeCollision.hh
LIB1_H += Collision/EdgeEdgeCollision.hh
LIB1_H += Collision/EdgeFaceCollision.hh
LIB1_H += Collision/EdgeFaceIntersection.hh
LIB1_H += Collision/ElementProxy.hh
LIB1_H += Collision/ImpulseBase.hh
LIB1_H += Collision/LinearCompliantImpulse.hh
LIB1_H += Collision/NonLinearCompliantImpulse.hh
LIB1_H += Collision/OrientedBoundingBox.hh
LIB1_H += Collision/ProximityCollision.hh
LIB1_H += Collision/TriangularMesh.hh
LIB1_H += Collision/VanillaImpulse.hh
LIB1_H += Collision/VertexFaceCollision.hh
LIB1_H += Core/BandMatrix.hh
LIB1_H += Core/BandMatrixFwd.hh
LIB1_H += Core/BoundingBox.hh
LIB1_H += Core/CollisionParameters.hh
LIB1_H += Core/Definitions.hh
LIB1_H += Core/ElasticStrand.hh
LIB1_H += Core/ElasticStrandParameters.hh
LIB1_H += Core/ElasticStrandUtils.hh
LIB1_H += Core/LinearSolver.hh
LIB1_H += Core/Resampler.hh
LIB1_H += Core/StepperBase.hh
LIB1_H += Core/StrandBackup.hh
LIB1_H += Core/StrandBase.hh
LIB1_H += Core/StrandState.hh
LIB1_H += Core/SurfaceInfo.hh
LIB1_H += Dependencies/ElasticityParameters.hh
LIB1_H += Dependencies/BendingProducts.hh
LIB1_H += Dependencies/DependencyNode.hh
LIB1_H += Dependencies/DegreesOfFreedom.hh
LIB1_H += Dependencies/Kappas.hh
LIB1_H += Dependencies/MaterialFrames.hh
LIB1_H += Dependencies/ReferenceFrames.hh
LIB1_H += Dependencies/Twists.hh
LIB1_H += Dynamic/ConstraintScriptingController.hh
LIB1_H += Dynamic/DOFScriptingController.hh
LIB1_H += Dynamic/ImplicitStepper.hh
LIB1_H += Dynamic/MeshScriptingController.hh
LIB1_H += Dynamic/ProblemFwd.hh
LIB1_H += Dynamic/ScriptingController.hh
LIB1_H += Dynamic/SimulationParameters.hh
LIB1_H += Dynamic/StrandImplicitManager.hh
LIB1_H += Dynamic/StrandDynamicTraits.hh
LIB1_H += Forces/AdhesiveAnchorForce.hh
LIB1_H += Forces/AirDragForce.hh
LIB1_H += Forces/AngularSpringForce.hh
LIB1_H += Forces/BendingForce.hh
LIB1_H += Forces/Bridson/Vector.hh
LIB1_H += Forces/Bridson/array1.hh
LIB1_H += Forces/Bridson/array3.hh
LIB1_H += Forces/Bridson/util.hh
LIB1_H += Forces/Bridson/vec.hh
LIB1_H += Forces/Clump.hh
LIB1_H += Forces/ClumpingForce.hh
LIB1_H += Forces/CombForce.hh
LIB1_H += Forces/ExplicitForce.hh
LIB1_H += Forces/ForceAccumulator.hh
LIB1_H += Forces/ForceBase.hh
LIB1_H += Forces/GaussianVolumetricForce.hh
LIB1_H += Forces/GravitationForce.hh
LIB1_H += Forces/InertialForce.hh
LIB1_H += Forces/LevelSet.hh
LIB1_H += Forces/LevelSetForce.hh
LIB1_H += Forces/LevelSetFwd.hh
LIB1_H += Forces/MassDampingForce.hh
LIB1_H += Forces/PressureForce.hh
LIB1_H += Forces/SpringForce.hh
LIB1_H += Forces/StrandStrandForce.hh
LIB1_H += Forces/StretchingForce.hh
LIB1_H += Forces/SurfaceConstraintForce.hh
LIB1_H += Forces/TwistingForce.hh
LIB1_H += Forces/ViscousOrNotViscous.hh
LIB1_H += Render/Color.hh
LIB1_H += Render/OpenGLDecl.hh
LIB1_H += Render/OpenGLHeaders.hh
LIB1_H += Render/RenderBase.hh
LIB1_H += Render/StrandRenderer.hh
LIB1_H += Static/CollisionSet.hh
LIB1_H += Static/CoupledStaticStepper.hh
LIB1_H += Static/GradientDescentStepper.hh
LIB1_H += Static/LogBucketQueue.hh
LIB1_H += Static/ReverseStaticStepper.hh
LIB1_H += Static/StaticCollision.hh
LIB1_H += Static/StaticCollisionsInfo.hh
LIB1_H += Static/StaticStepper.hh
LIB1_H += Static/StaticForceList.hh
LIB1_H += Static/StrandStaticManager.hh
LIB1_H += Static/StrandStaticTraits.hh
LIB1_H += Static/StrandStrandForcesManager.hh
LIB1_H += Static/SynchronizationPoint.hh
LIB1_H += Utils/Cycle.hh
LIB1_H += Utils/Distances.hh
LIB1_H += Utils/EigenSerialization.hh
LIB1_H += Utils/IndexedSet.hh
LIB1_H += Utils/LoggingTimer.hh
LIB1_H += Utils/Memory.hh
LIB1_H += Utils/MetricTree.hh
LIB1_H += Utils/MsgInfo.hh
LIB1_H += Utils/SpatialHashMap.hh
LIB1_H += Utils/SpatialHashMapFwd.hh
LIB1_H += Utils/SymmetricBandMatrixSolver.hh
LIB1_H += Utils/TextLog.hh
LIB1_H += Utils/ThreadUtils.hh
LIB1_H += Utils/TimeUtils.hh
LIB1_H += Utils/Timer.hh

# Source files
LIB1_SDIR := .
LIB1_S += Collision/BVH.cc
LIB1_S += Collision/CollisionBase.cc
LIB1_S += Collision/CollisionDetector.cc
LIB1_S += Collision/CollisionUtils.cc
LIB1_S += Collision/ContinuousTimeCollision.cc
LIB1_S += Collision/EdgeEdgeCollision.cc
LIB1_S += Collision/EdgeFaceCollision.cc
LIB1_S += Collision/EdgeFaceIntersection.cc
LIB1_S += Collision/ElementProxy.cc
LIB1_S += Collision/ImpulseBase.cc
LIB1_S += Collision/LinearCompliantImpulse.cc
LIB1_S += Collision/NonLinearCompliantImpulse.cc
LIB1_S += Collision/OrientedBoundingBox.cc
LIB1_S += Collision/ProximityCollision.cc
LIB1_S += Collision/TriangularMesh.cc
LIB1_S += Collision/VanillaImpulse.cc
LIB1_S += Collision/VertexFaceCollision.cc
LIB1_S += Core/ElasticStrand.cc
LIB1_S += Core/ElasticStrandParameters.cc
LIB1_S += Core/ElasticStrandUtils.cc
LIB1_S += Core/LinearSolver.cc
LIB1_S += Core/StrandBackup.cc
LIB1_S += Core/StrandBase.cc
LIB1_S += Core/StrandState.cc
LIB1_S += Core/SurfaceInfo.cc
LIB1_S += Dependencies/BendingProducts.cc
LIB1_S += Dependencies/DegreesOfFreedom.cc
LIB1_S += Dependencies/Kappas.cc
LIB1_S += Dependencies/MaterialFrames.cc
LIB1_S += Dependencies/ReferenceFrames.cc
LIB1_S += Dependencies/Twists.cc
LIB1_S += Dynamic/ConstraintScriptingController.cc
LIB1_S += Dynamic/DOFScriptingController.cc
LIB1_S += Dynamic/ImplicitStepper.cc
LIB1_S += Dynamic/MeshScriptingController.cc
LIB1_S += Dynamic/ScriptingController.cc
LIB1_S += Dynamic/SimulationParameters.cc
LIB1_S += Dynamic/StrandImplicitManager.cc
LIB1_S += Dynamic/StrandDynamicTraits.cc
LIB1_S += Forces/AdhesiveAnchorForce.cc
LIB1_S += Forces/AirDragForce.cc
LIB1_S += Forces/AngularSpringForce.cc
LIB1_S += Forces/BendingForce.cc
LIB1_S += Forces/Clump.cc
LIB1_S += Forces/ClumpingForce.cc
LIB1_S += Forces/CombForce.cc
LIB1_S += Forces/ExplicitForce.cc
LIB1_S += Forces/ForceBase.cc
LIB1_S += Forces/GaussianVolumetricForce.cc
LIB1_S += Forces/GravitationForce.cc
LIB1_S += Forces/InertialForce.cc
LIB1_S += Forces/LevelSet.cc
LIB1_S += Forces/LevelSetForce.cc
LIB1_S += Forces/MassDampingForce.cc
LIB1_S += Forces/PressureForce.cc
LIB1_S += Forces/SpringForce.cc
LIB1_S += Forces/StrandStrandForce.cc
LIB1_S += Forces/StretchingForce.cc
LIB1_S += Forces/SurfaceConstraintForce.cc
LIB1_S += Forces/TwistingForce.cc
LIB1_S += Render/RenderBase.cc
LIB1_S += Render/StrandRenderer.cc
LIB1_S += Static/CoupledStaticStepper.cc
LIB1_S += Static/CollisionSet.cc
LIB1_S += Static/GradientDescentStepper.cc
LIB1_S += Static/LogBucketQueue.cc
LIB1_S += Static/ReverseStaticStepper.cc
LIB1_S += Static/StaticCollisionsInfo.cc
LIB1_S += Static/StaticStepper.cc
LIB1_S += Static/StrandStaticManager.cc
LIB1_S += Static/StrandStaticTraits.cc
LIB1_S += Static/StrandStrandForcesManager.cc
LIB1_S += Utils/Distances.cc
LIB1_S += Utils/Memory.cc
LIB1_S += Utils/MsgInfo.cc
LIB1_S += Utils/TextLog.cc
LIB1_S += Utils/Timer.cc


##############################################################################
include /vol/bob/make/tail.mk
##############################################################################
