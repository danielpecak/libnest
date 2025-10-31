# 📊 ANALIZA PROJEKTU LIBNEST

**Typ projektu**: Biblioteka naukowa Python do fizyki jądrowej i gwiazd neutronowych  
**Rozmiar**: ~5000 linii kodu Python  
**Status**: Projekt badawczy z dokumentacją Sphinx  
**Data analizy**: 1 listopada 2025

---

## 🎯 PLAN POPRAWEK I ULEPSZEŃ

### **1. KRYTYCZNE - Struktura Projektu**

#### ❌ **Problem**: Brak pliku konfiguracyjnego pakietu
- Brak `setup.py` lub `pyproject.toml`
- Niemożliwa instalacja przez `pip install`
- Brak metadanych pakietu

#### ✅ **Rozwiązanie**:
- Utworzyć `pyproject.toml` z konfiguracją Poetry/setuptools
- Dodać informacje o licencji, autorze, zależnościach
- Umożliwić instalację: `pip install -e .`

---

### **2. WYSOKIEJ WAGI - Jakość Kodu**

#### ❌ **Problemy**:
1. **`main.py`** - zaśmiecony plik testowy (6674 bajty zakomentowanego kodu)
2. **Brak testów jednostkowych** - w `tests/` są tylko podstawowe testy
3. **TODOs w kodzie** - 21 nierozwiązanych zadań
4. **Słaba dokumentacja README.md** - tylko "# libnest"
5. **Niekonsystentne nagłówki** - mieszanka polskich i angielskich komentarzy

#### ✅ **Rozwiązania**:
- Wyczyścić `main.py` lub przenieść do `examples/`
- Napisać właściwe testy jednostkowe z pytest
- Rozwiązać lub udokumentować wszystkie TODOs
- Napisać kompletny README z przykładami użycia
- Ujednolicić język (preferowane EN dla kodu)

---

### **3. ŚREDNIEJ WAGI - Jakość Kodu**

#### ❌ **Problemy w `libnest/tools.py`**:
```python
# Linia 134-135 - niezdefiniowane zmienne globalne
energy=e*HBARC/(2*(MN*particleN(density_n)+MP*particleN(density_p))/(particleN(density)))
```
- `HBARC`, `MN`, `MP` nie są zaimportowane

#### ❌ **Problem w `libnest/definitions.py`**:
- Linia 105, 301, 339 - TODOs o połączeniu funkcji
- Funkcja `mu_q()` używa `effMn()` i `effMp()` które nie istnieją w pliku

#### ✅ **Rozwiązania**:
- Dodać brakujące importy z `libnest.units`
- Połączyć duplikujące się funkcje
- Sprawdzić i naprawić wszystkie zależności

---

### **4. ŚREDNIEJ WAGI - Dokumentacja**

#### ❌ **Problemy**:
- Minimalistyczny README
- Brak przykładów użycia w głównym katalogu
- Brak informacji o wymaganiach systemowych
- Brak badge'ów (testy, pokrycie, wersja)

#### ✅ **Rozwiązania**:
- Rozbudować README o:
  - Opis projektu i celów naukowych
  - Instrukcję instalacji
  - Szybki start (quick start)
  - Przykłady użycia
  - Link do pełnej dokumentacji
  - Informacje o autorach i finansowaniu
  - Status rozwoju

---

### **5. NISKIEJ WAGI - Best Practices**

#### ❌ **Problemy**:
- Brak CI/CD (GitHub Actions)
- Brak formattera (black, ruff)
- Brak pre-commit hooks
- Brak type hints
- Kodowanie `# -*- coding:utf-8 -*-` niepotrzebne w Python 3

#### ✅ **Rozwiązania**:
- Dodać GitHub Actions dla testów
- Skonfigurować `ruff` lub `black`
- Dodać pre-commit hooks
- Stopniowo dodawać type hints
- Usunąć przestarzałe nagłówki

---

### **6. OPTYMALIZACJE**

#### Możliwe usprawnienia w `libnest/tools.py`:
- Funkcja `threeSlice()` używa `dtype="object"` - lepiej używać konkretnych typów
- `centerOfMass()` - można zoptymalizować używając `np.average()` z wagami

---

## 📋 PRIORYTETOWA LISTA ZADAŃ

### **FAZA 1: Fundamenty** (1-2 dni)
1. ✅ Utworzenie `pyproject.toml`
2. ✅ Napisanie kompletnego README.md
3. ✅ Dodanie `setup.py` lub pełna konfiguracja pyproject.toml
4. ✅ Porządkowanie `main.py` (przeniesienie do examples/)

### **FAZA 2: Jakość** (2-3 dni)
5. ✅ Naprawienie brakujących importów
6. ✅ Rozwiązanie TODOs w kodzie
7. ✅ Napisanie testów jednostkowych
8. ✅ Dodanie GitHub Actions CI

### **FAZA 3: Polerowanie** (1-2 dni)
9. ✅ Konfiguracja formattera (ruff/black)
10. ✅ Dodatnie type hints
11. ✅ Usprawnienie dokumentacji docstring
12. ✅ Dodanie przykładów użycia

---

## 📝 SZCZEGÓŁOWE ZNALEZIONE PROBLEMY

### TODOs w kodzie (21 wystąpień):

#### `libnest/myio.py` (linia 19)
- Handling WDATA format with Forbes' python library

#### `libnest/bsk.py`
- Linia 18: Make table of BSk parameters
- Linia 567, 732: TODO lists w komentarzach
- Linia 999: Make different functions for neutrons and protons

#### `libnest/nucleus.py`
- Linia 50, 123: DEALING WITH background density
- Linia 62: Include X, Y, Z directions
- Linia 98: DOCUMENTATION, FORMULA

#### `libnest/definitions.py`
- Linia 105: Check if I should return rho_n-rho_p OR rho_n-rho_p/(rho_n+rho_p + DENSEPSILON)
- Linia 301, 339: Join these two functions: E_minigap_delta_n and E_minigap_rho_n

#### `main.py`
- Linia 9: TODO - general file description

---

## 🔍 BRAKUJĄCE ZALEŻNOŚCI

### `libnest/tools.py`
- Linia 134-135: Używa `HBARC`, `MN`, `MP` bez importu z `libnest.units`

### `libnest/definitions.py`
- Funkcja `mu_q()` wywołuje `effMn()` i `effMp()` - prawdopodobnie z `libnest.bsk`
- Brak importu `sys` mimo używania `sys.exit()`

---

## 📦 STRUKTURA DEPENDENCIES

### Obecne (requirements.txt):
- sphinxcontrib-bibtex
- sphinx_rtd_theme
- matplotlib
- numpy
- pandas

### Sugerowane dodatkowe:
- pytest (testy)
- pytest-cov (pokrycie testów)
- ruff lub black (formatowanie)
- mypy (type checking)
- pre-commit (hooks)

---

## 🎓 KONTEKST NAUKOWY

Projekt jest częścią prac badawczych z zakresu:
- Fizyki jądra atomowego
- Struktury gwiazd neutronowych
- Parametryzacji Brussels-Montreal (BSk)
- Teorii parowania i nadciekłości neutronów

**Afiliacje**:
- Warsaw Technical University
- Université Libre de Bruxelles
- Institute of Physics, Polish Academy of Sciences

**Autor**: Daniel Pęcak <daniel.pecak@pw.edu.pl>

---

## 📊 STATYSTYKI KODU

- **Całkowita liczba linii**: ~5011
- **Liczba modułów**: 8 głównych modułów w `libnest/`
- **Liczba skryptów przykładowych**: 1 (`examples/tools_example.py`)
- **Liczba skryptów do dokumentacji**: 11 (w `docs/source/plots/`)
- **Liczba plików testowych**: 2 (`tests/tests.py`, `tests/units_tests.py`)

---

## ✅ NASTĘPNE KROKI

1. **Natychmiastowe**: Naprawienie brakujących importów
2. **Krótkoterminowe**: Czyszczenie kodu, rozwiązanie TODOs
3. **Średnioterminowe**: Dodanie testów i CI/CD
4. **Długoterminowe**: Dokumentacja i publikacja pakietu
