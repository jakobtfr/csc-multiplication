# Projektbericht

## Aufgabenstellung

In diesem Projekt wurde die Multiplikation dünnbesetzer Matrizen mit dem
Compressed Sparse Column (CSC)-Format implementiert.
Hier werden die Dimensionen der Matrix, die nicht-null Werte (values), die Zeilenindizes (row_ind) 
und die Spaltenpointer (colPtr), 
die anzeigen, wo eine neue Spalte beginnt, gespeichert, um Speichereffizienz zu erhöhen.

## Lösungsansätze

### Naiver Ansatz (V2)
Die traditionelle Matrixmultiplikation, bei der das Skalarprodukt zwischen den Zeilen von A und den Spalten von B berechnet wird. 
Dies erfordert ineffizientes Iterieren über jede Spalte von A.
### Optimierung 1 (V1)
Effizientere Multiplikation durch Transponieren der Matrix: Berechnung des Skalarprodukts der Spalten von A mit den Spalten von B.
### Optimierung 2 (V2)
Reduzierung unnötiger Operationen durch in-place Optimierung, vermeidet das Finden und Kopieren der Spalten in V1 von A und B.
## Benchmarking
### Methodik
Zufällig generierte Matrizen mit Dichten zwischen 10% und 33% wurden verwendet. Einmalig erstellte Matrizen wurden für alle Tests verwendet.
- Laufzeit: Es wurden die Produkte von 32x32, 64x64, ..., 8192x8192 Matrizen berechnet. Für kleinere
Matrizen (32x32 bis 512x512) wurde öfter iteriert und der Durchschnitt berechnet. Größere Matrizen
einmalig.
  - 32x32: 10.000 Iterationen
  - 64x64: 10.000 Iterationen
  - 128x128: 1.000 Iterationen
  - 256x256: 100 Iterationen
  - 512x512: 100 Iterationen
- Profiling: Perf für das Produkt zweier 2048x2048 Matrizen
- Heap-Profiling: Massif für das Produkt zweier 1024x1024
Matrizen, da Laufzeiterhöhung durch Massif\

### Testumgebung
- CPU: Core i9-13900H
- Speicher: 16 GB 
- Betriebssystem: Fedora Linux Workstation 40\
- Compiler: GCC 14.1.1 
- CFLAGS: -O2, -Wall -Wextra -Wpedantic

Profiling CFLAGS: -g -fno-omit-frame-pointer\

Heap-Profiling CFLAGS: -g

### Transponieren der Matrix
Transponierung mit RadixSort\
Mit Sortieren der row_Ind (col_Ind der Ursprungsmatrix) kann A transponiert werden.
Gleichzeitige Sortierung von values und col_ind.

### Parsing und Rahmenprogramm
Rahmenprogramm\
Das Rahmenprogramm ermöglicht die Bestimmung verschiedener Optionen:
- Argumente der Eingabe-Matrizen, Ausgabedatei
- Versionsauswahl
- Zeitmessung
- Programmwiederholungen
- Echtzeit-Überwachung des Programmfortschritts

Parsing\
Gleichzeitige Erstellung von nötigen Daten für transpose und Verarbeitung der Matrizen.

### Multiplikation
Es wurden drei Versionen des Multiplikationsalgorithmus implementiert. Gemeinsame Entscheidung, 
die Multiplikation der transponierten A-Matrix und der B-Matrix
zu implementieren, da diese wahrscheinlich performanter ist.\
Die originale Implementierung war nicht optimal, wegen unnötigem Allokieren und Kopieren von
Werten.
In V2 alles möglichst In-Place, um dies zu vermeiden.\
Zuletzt wurde der naive Ansatz implementiert. Damit
konnte sicherstellt werden, dass die ersten zwei Implementierungen tatsächlich performanter sind.\

Sonstiges: Implementierung von Makefile, Schreiben von Tests und Hilfsfunktionen\

