(* ::Package:: *)

(* ::Section:: *)
(*Beginning*)


BeginPackage["QSWalk`"];


QuantumStochasticWalk::usage = "QuantumStochasticWalk[H, LkList, omega, rho0, t] computes a quantum stochastic walk \
with Hamiltonian matrix H, Lindblad operators LkList, and weighting factor omega; starting from density matrix rho0 at \
time 0 and returning the density matrix at time t. H should be an NxN Hermitian matrix, and LkList should be a list of \
NxN matrices; rho0 should be an NxN Hermitian matrix with a trace of unity.";


QuantumWalk::usage = "QuantumWalk[H, psi0, t] computes a (continuous-time) quantum walk \
with Hamiltonian matrix H, starting from state vector psi0 at time 0 and returning the \
state vector at time t. H should be an NxN Hermitian matrix; psi0 should be a length-N list of \
complex-valued probability amplitudes with squared magnitudes summing to unity.";


ClassicalRandomWalk::usage = "ClassicalRandomWalk[M, p0, t] computes a (continuous-time) classical random walk \
with generator matrix M, starting from probability vector p0 at time 0 and returning the probability vector at time t. \
M should be a real-valued NxN matrix; p0 should be a length-N list of non-negative probabilities summing to unity.";


GeneratorMatrix::usage = "GeneratorMatrix[G, gamma] returns the generator matrix M with rate parameter gamma \
for the graph G.";


LindbladSet::usage = "LindbladSet[mat] returns of set of Lindblad matrices corresponding to each non-zero element of mat.";


GoogleMatrix::usage = "GoogleMatrix[G, alpha] returns the \"Google matrix\" with damping factor alpha for the graph G.";


Vectorize::usage = "Vectorize[mat] returns a vector formed by concatenating the columns of matrix mat.";


Matricize::usage = "Matricize[vec] returns a matrix formed by inverting the action of the Vectorize function. \
The length of vec must be an integer squared.";


CayleyTree::usage = "CayleyTree[d, n] returns a Graph object representing an n-th generation Cayley tree of order d.";


GluedBinaryTree::usage = "GluedBinaryTree[n] returns a graph comprising two complete, (n+1)-level binary trees glued together \
in order along their leaf nodes.";


RandomGluedBinaryTree::usage = "RandomGluedBinaryTree[n] returns a graph comprising two complete, (n+1)-level binary trees \
glued together in random order along their leaf nodes.";


Begin["`Private`"];


(* ::Section:: *)
(*Utility functions*)


Vectorize[mat_?MatrixQ] := Flatten@Transpose[mat]


Matricize[vec_?VectorQ] := Transpose@Partition[vec, Sqrt[Length[vec]]] /; IntegerQ[Sqrt[Length[vec]]]


LindbladSet[mat_?MatrixQ] := SparseArray[{#}, Dimensions[mat]] & /@ Most[ArrayRules[Sqrt[Abs[mat]]]]


GeneratorMatrix[G_Graph, gamma_] := Module[{A=Transpose[WeightedAdjacencyMatrix[G]]},
	gamma (DiagonalMatrix[SparseArray@Total[A]] - A)]


GoogleMatrix[A_?MatrixQ, alpha_] := Module[{n=Length[A], tot=Total[A], one=N[1,Precision[A]], fill, dnom},
	fill = one - Unitize[tot];
	dnom = tot + n*fill;
	alpha Map[(#+fill)/dnom &, A] + (one-alpha)/n
]


GoogleMatrix[G_Graph, alpha_] := GoogleMatrix[Transpose@WeightedAdjacencyMatrix[G], alpha]


(* ::Section:: *)
(*Graphs*)


Options[CayleyTree]={DirectedEdges->False};
CayleyTree[k_Integer /; k>1, level_Integer?Positive, opts:OptionsPattern[]] := Module[{edgeFn, edgeList},
edgeFn = If[OptionValue[DirectedEdges], DirectedEdge, UndirectedEdge];
edgeList = Flatten @ NestList[Table[edgeFn[#[[1+Quotient[j-1,k-1],2]], #[[-1,2]]+j], {j,(k-1)*Length[#]}] &,
	Table[edgeFn[1,j+1], {j,k}], level-1];
Graph[edgeList, GraphLayout->"RadialDrawing", Sequence@@FilterRules[{opts}, Options[Graph]]]]


GluedBinaryTree[n_, opts:OptionsPattern[]] := Module[{M=2^n+2^(n+1)-2, leftEdges, rightEdges},
	(* REF: http://arxiv.org/abs/quant-ph/0209131v2, Fig. 1 *)
	leftEdges = EdgeList @ CompleteKaryTree[n+1];
	rightEdges = leftEdges /. UndirectedEdge[i_,j_] :> UndirectedEdge[M-j+1,M-i+1];
	Graph[Range[M], Join[leftEdges, rightEdges], Sequence@@FilterRules[{opts}, Options[Graph]]]]


RandomGluedBinaryTree[n_, opts:OptionsPattern[]] := Module[{M=2^(n+2)-2, leftEdges, rightEdges, leftLeaves, rightLeaves, centreCycle, centreEdges},
	(* REF: http://arxiv.org/abs/quant-ph/0209131v2, Fig. 2 *)
	leftEdges = EdgeList @ CompleteKaryTree[n+1];
	rightEdges = Reverse[leftEdges] /. UndirectedEdge[i_,j_] :> UndirectedEdge[M-j+1,M-i+1];
	leftLeaves = Range[(M+2)/4, M/2];
	rightLeaves = Range[M/2+1, (3M+2)/4];
	centreCycle = Riffle[RandomSample[leftLeaves], RandomSample[rightLeaves]];
	centreEdges = UndirectedEdge @@@ Partition[centreCycle, 2, 1, {1,1}];
	Graph[Range[M], Join[leftEdges, centreEdges, rightEdges], Sequence@@FilterRules[{opts}, Options[Graph]]]]


(* ::Section:: *)
(*Walk functions*)


(* ::Subsection:: *)
(*Classical random walk*)


ClassicalRandomWalk[M_, p0_, t_] := MatrixExp[-M t, p0]


(* ::Subsection:: *)
(*Quantum walk*)


QuantumWalk[H_, psi0_, t_] := MatrixExp[-I H t, psi0]


(* ::Subsection:: *)
(*Quantum stochastic walk*)


QuantumStochasticWalk[H_, LkList_, omega_, rho0_, t_] := Module[{n=Length[H], id, LdagL, masterL},
	id = IdentityMatrix[n, SparseArray];
	(* vectorized master equation matrix *)
	masterL = -(1-omega) I (KroneckerProduct[id, H] - KroneckerProduct[Transpose[H],id]) +
		omega Sum[KroneckerProduct[Conjugate[L], L] - (1/2)(KroneckerProduct[id, ConjugateTranspose[L].L] + 
			KroneckerProduct[Transpose[L].Conjugate[L], id]), {L, LkList}];
	(* solution *)
	Matricize @ MatrixExp[masterL t, Vectorize[rho0]]
]


(* ::Section:: *)
(*Ending*)


End[];


EndPackage[];
