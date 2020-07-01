Needs["QSWalk`"]

NaturalSort[list_]:= list[[Ordering @ PadRight @ StringSplit[list, x : DigitCharacter .. :> FromDigits @ x]]]

SO[H_, LkList_, omega_] := Module[{n=Length[H], id, LdagL, masterL},
	id = IdentityMatrix[n, SparseArray];
	(* vectorized master equation matrix *)
	masterL = -(1-omega) I (KroneckerProduct[id, H] - KroneckerProduct[Transpose[H],id]) +
		omega Sum[KroneckerProduct[Conjugate[L], L] - (1/2)(KroneckerProduct[id, ConjugateTranspose[L].L] +
			KroneckerProduct[Transpose[L].Conjugate[L], id]), {L, LkList}];
	Return[masterL]
]

stepbench[omega_, time_, files_, filessym_, results_, maxsimtime_]:=
	Module[{totalsimtime, G, Gsym, totalsteptime, H, L, rho0, rhot, file},
		totalsimtime = 0.0;
		For[i = 1, i < Length[files], i++,

			G = Import[files[[i]]];
			Gsym = Import[filessym[[i]]];

			totalsteptime = SessionTime[];

			H = Gsym;

			L = LindbladSet[G^2];

            rho = N[IdentityMatrix[Length[G]]/Length[G]];

			sotime = SessionTime[];
			masterL = SO[Gsym, L, omega];
			sotime = SessionTime[] - sotime;

			exptime = SessionTime[];
			rhot = Matricize[MatrixExp[masterL t, Vectorize[rho]]];
			exptime = SessionTime[] - exptime;

			totalsteptime  = SessionTime[] - totalsteptime;
			totalsimtime = totalsimtime + totalsteptime;

			peakmemory = N[MaxMemoryUsed[]/1073741824];

			savetime = SessionTime[];
			Export[results, FileBaseName[files[[i]]] -> rhot, OverwriteTarget -> "Append"];
			savetime = SessionTime[] - savetime;

			file = OpenAppend["Results/QSWalk_m_local_step.csv"];
			WriteString[file, StringForm["`1`,`2`,`3`,`4`,`5`,`6`,`7`,`8`,`9`\n",
			FileBaseName[files[[i]]],
			Length[G],
			Length[G["NonzeroPositions"]],
			Length[masterL["NonzeroPositions"]],
			sotime,
			exptime,
			totalsteptime,
			peakmemory,
			savetime]];
			Close[file];

			If[totalsimtime > maxsimtime, Break[]]
		]
	]

(* Overwrite durring first run (lines)*)

If[! FileExistsQ["Results/QSWalk_m_local_step.csv"],
	timesfile = OpenWrite["Results/QSWalk_m_local_step.csv"];
	WriteString[timesfile, "name,dim,nnz,SO_nnz,SO_time,step_time,total_time,peak_memory,save_time\n"];
	Close[timesfile];
]


simtime = ToExpression[$ScriptCommandLine[[5]]]
omega = 0.1
t = 100.0

files = NaturalSort[FileNames[$ScriptCommandLine[[2]] <> "*.mtx"]]
filessym = NaturalSort[FileNames[$ScriptCommandLine[[3]] <> "*.mtx"]]
stepbench[omega, t, files, filessym, $ScriptCommandLine[[4]] <> "_QSWalk_m.h5", simtime]
