BeginPackage["Quantum`"]
(*Aquí va la documentación*)
GeneradorQs::usage = "GeneradorQs recibe la lista de proyecciones
 de espín del estado y devuelve los valores
 de proyección intermedios:
 GeneradorQs : {Semienteros} -> {Semienteros}
 Ejemplo : {1/2, 1/2, 1/2} -> {1}."

RangoAngular::usage = "RangoAngular devuelve los valores de proyección
para un valor de espín dado.
RangoAngular: Semientero -> {Semienteros}
Ejemplo: 1 -> {-1,0,1}"

QVD::usage = " QVD nos dice si los valores de Q son válidos, i.e., 
corresponden a valores válidos de la proyección de los
momentos angulares intermedios.
QVD: {Semienteros}, {Semienteros} -> Bool
Ejemplos: {1/2,1}, {1/2,1/2} -> False"

Mexico::usage = "México toma los valores de espín intermedios, sus proyecciones,
proyecciones de los estados, j y m. Devuelve el valor
del coeficiente de Clebsch-Gordan generalizado.
México: {Semientero},{Semientero},{Semientero},Semientero,Semientero -> Número real"

Lapiz::usage = " Lapiz toma los valores de momento angular intermedios,
las proyecciones de los estados, j y m. Devuelve dos
valores booleanos que nos dicen si las proyecciones
son consistentes con el valor de m que se le dió y
si las proyecciones que se calculan con GeneradorQs
son válidas con los valores de espín intermedios.
Lapiz: {Semienteros}, {Semienteros}, Semientero, Semientero ->
{Semienteros}"

Estados::usage = "Estados genera la base computacional a partir del
número de qbits.
Estados: Entero -> {Semienteros}
Ejemplo: Estados[2]: {1/2,1/2}, {1/2, -1/2}, {-1/2,1/2}, {-1/2, -1/2}
"

PasoEspinAEstado::usage = "PasoEspinAEstado está pensado para ser
utilizado únicamente con Map.
Map le pasa a PasoEspinAEstado cada una de las entradas de una lista
con 1/2 o -1/2 y la convierte en {1,0} o {0,1} para generar los
eigenestados de la matriz de Pauli, z.
Lapiz: Semientero(1/2 o -1/2) -> {0,1} o {1,0}
"
EspinAEstado::usage = "Genera el estado utilizando PasoEspinAEstado.
EspinAEstado: {1/2 o -1/2} -> Estado base computacional."
Begin["Private`"]
(*Aquí van las funciones*)

GeneradorQs[M_] := Module[{l = Length[M], q, i},
  q = {M[[l]] + M[[l - 1]]};
  Table[AppendTo[q, q[[i]] + M[[l - (i + 1)]]], {i, 1, l - 3}];
  q]

RangoAngular[j_] := Table[m, {m, -j, j}]

QVD[k_, q_] := Module[{i, total = 0},
  For[i = 1, i <= Length[k], i++, 
   If[MemberQ[RangoAngular[k[[i]]], q[[i]]], total += 1]] ; 
  If[total == Length[k], True, False]]

Mexico[k_, q_, M_, j_, m_] := Module[{l = Length[M], prod = 0, i},
  prod = ClebschGordan[{1/2, M[[l - 1]]}, {1/2, M[[l]]}, {k[[1]], 
     q[[1]]}];
  For[i = 1, i <= l - 3, i++, 
   prod *= ClebschGordan[{1/2, M[[l - i - 1]]}, {k[[i]], 
      q[[i]]}, {k[[i + 1]], q[[i + 1]]}]];
  prod *= 
   ClebschGordan[{1/2, M[[1]]}, {k[[l - 2]], q[[l - 2]]}, {j, m}]]

Lapiz[K_, M_, j_, m_] := If[Total[M] == m && QVD[K, GeneradorQs[M]],
  Mexico[K, GeneradorQs[M], M, j, m]*EspinAEstado[M],
  Table[0,{i,2^Length[M]}]]

Estados[n_]:=Tuples[{1/2,-1/2},{n}]

PasoEspinAEstado[estado_] := If[estado == 1/2, {1, 0}, {0, 1}]

EspinAEstado[estados_] := 
 Flatten[KroneckerProduct @@ Map[PasoEspinAEstado, estados]]
End[] 
EndPackage[]
