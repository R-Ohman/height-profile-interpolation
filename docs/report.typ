#set heading(
  numbering: "1.",
)

#show heading.where(level: 1): set block(above: 1.75em, below: 1.25em)
#show heading.where(level: 1): set text(size: 14pt)

#show figure.caption: it => [
  #if it.supplement == [Figure] [
    Rysunek  #it.counter.display()#it.separator #it.body
  ] else if it.supplement == [Table] [
    Tabela  #it.counter.display()#it.separator #it.body
  ] else [
    #it.supplement #it.numbering#it.separator #it.body
  ]
]

#set page(
  paper: "a4",
  margin: 2cm,
  footer: context [
    #set text(8pt)
    *Projekt 3* - sprawozdanie
    #h(1fr)
    #if counter(page).get().first() > 1 [
      #counter(page).display(
        "1 / 1",
        both: true
      )
    ]
  ]
)

#let creationDate = datetime.today()
#align(right)[
    #stack(
      image("plots/pg.jpg", width: 15%)
    )
  ]
#align(center)[
  #stack(
    dir: ttb,
    text(size: 24pt, weight: "semibold")[Aproksymacja profilu wysokościowego],
    v(10pt),
    text(size: 14pt)[Metody Numeryczne - Projekt 3],
    v(20pt),
    text(size: 12pt)[Ruslan Rabadanov 196634],
    v(10pt),
    text(size: 12pt)[#creationDate.display("[day]-[month]-[year]")]
  )
]

#v(40pt)

#outline(
  title: "Spis treści",
)

#v(30pt)

= Wstęp

Celem tego projektu jest implementacja i porównanie algorytmów interpolacji funkcji, używając przykładu aproksymacji profilu wysokościowego terenu.

W projekcie zaimplementowano metodę interpolacji wielomianowej Lagrange'a i metodę interpolacji z użyciem funkcji sklejanych trzeciego stopnia
wraz z niezbędnymi funkcjami macierzowymi i algorytmami do rozwiązywania układów równań liniowych.

Implementacja została przeprowadzona w języku `Python`, bez korzystania z zewnętrznych bibliotek do obliczeń numerycznych. Wykorzystano bibliotekę `matplotlib` do rysowania wykresów.

Projekt ma na celu porównanie wyników interpolacji dla różnych zestawów danych, charakteryzujących się różnymi cechami, aby zbadać zachowanie obu metod w różnych warunkach. Dodatkowo, przeanalizowane zostaną różnice wynikające z różnych parametrów interpolacji, takich jak charakter funkcji, liczba i rozmieszczenie węzłów interpolacji.

== Interpolacja wielomianowa Lagrange'a

Interpolacja wielomianowa Lagrange'a polega na znalezieniu wielomianu stopnia n-1, który przechodzi przez n punktów. 

W tym celu przyda się policzenie bazy Lagrange'a:

$ phi.alt_i (x) = limits(product)_(j=1,j!=i)^(n+1) frac(x - x_j,x_i - x_j) , i = 1, 2, ..., n+1 $

Funkcja interpolująca ma postać:

$ F(x) = sum_(i=1)^(n+1) y_i phi.alt_i (x) $

Metoda ta jest stabilna i łatwa do implementacji, jednak ma swoje wady w postaci efektu Rungego, który polega na oscylacjach wynikowego wielomianu na krańcach przedziału interpolacji.

Efekt ten może zostać zniwelowany (jednak nie całkowicie wyeliminowany) poprzez zmianę rozmieszczenia węzłów interpolacji. W projekcie zaimplementowano dwa rodzaje węzłów: równoodległe oraz Czebyszewa.

== Interpolacja funkcjami sklejanymi trzeciego stopnia

Interpolacja funkcjami sklejanymi trzeciego stopnia polega na znalezieniu funkcji sklejanej, która przechodzi przez wszystkie węzły, a jej pochodne pierwszego i drugiego rzędu są ciągłe.

W _i_-tym podprzedziale wielomian sklejanej funkcji trzeciego stopnia ma postać:

$ S_i (x) = a_i + b_i (x - x_i) + c_i (x - x_i)^2 + d_i (x - x_i)^3 $

gdzie współczynniki $a_i, b_i, c_i, d_i$ są znane dla każdego z przedziałów.

Metoda ta jest bardziej skomplikowana w implementacji, jednak rozwiazuje problem efektu Rungego, który występuje w przypadku interpolacji wielomianowej Lagrange'a. Funkcje sklejane pozwalają na uzyskanie gładkiej interpolacji, która nie wykazuje oscylacji na krańcach przedziału interpolacji.

== Zestawy danych

Do analizy wybrano trzy zestawy danych, które różnią się swoimi cechami, aby zbadać zachowanie obu metod interpolacji w różnych warunkach:
- `chelm.txt` - płynna trasa z niewielkimi zmianami wysokości,
- `GlebiaChallengera.csv` - płaska trasa z jednym wzniesieniem,
- `100.csv` - trasa z wieloma wzniesieniami i spadkami.

#figure(
  image("plots/chelm/original.png"),
  caption: [Funkcja terenu dla zestawu danych Chelm]
)
#figure(
  image("plots/GlebiaChallengera/original.png"),
  caption: [Funkcja terenu dla zestawu danych Głębia Challengera]
)
#figure(
  image("plots/100/original.png"),
  caption: [Funkcja terenu dla zestawu danych 100]
)
// #pagebreak()

= Analiza podstawowa interpolacji

Analiza wyników interpolacji dla każdego zestawu danych koncentruje się na badaniu wpływu liczby punktów pomiarowych oraz ich rozmieszczenia na jakość interpolacji. Badanie zostało przeprowadzone dla zestawów danych `chelm.txt` oraz `GlebiaChallengera.csv`, uwzględniając od kilku do stu punktów pomiarowych spośród około pięciuset dostępnych.

Jak przedstawiono na poniższych wykresach, liczba punktów pomiarowych znacząco wpływa na jakość interpolacji. W przypadku interpolacji wielomianowej Lagrange'a z równoodległymi węzłami, efekt Rungego jest bardzo widoczny już przy 13 punktach pomiarowych, co znacznie ogranicza efektywność tej metody. Zastosowanie węzłów Czebyszewa pomaga uniknąć efektu Rungego do pewnej (dużej) liczby węzłów interpolacji. Przy użyciu funkcji sklejanych trzeciego stopnia, efekt Rungego nie występuje, a interpolacja jest dokładna nawet przy niewielkiej liczbie punktów pomiarowych.
#pagebreak()
== Interpolacja wielomianowa Lagrange'a

#figure(
  image("plots/chelm/linspace_nodes.png"),
  caption: [Wpływ liczby punktów pomiarowych na dokładność interpolacji wielomianowej dla zestawu danych Chelm (węzły równoodległe)]
)

#figure(
  image("plots/chelm/chebyshev_nodes.png"),
  caption: [Wpływ liczby punktów pomiarowych na dokładność interpolacji wielomianowej dla zestawu danych Chelm (węzły Czebyszewa)]
)

#figure(
  image("plots/GlebiaChallengera/linspace_nodes.png"),
  caption: [Wpływ liczby punktów pomiarowych na dokładność interpolacji wielomianowej dla zestawu danych Głębia Challengera (węzły równoodległe)]
)

#figure(
  image("plots/GlebiaChallengera/chebyshev_nodes.png"),
  caption: [Wpływ liczby punktów pomiarowych na dokładność interpolacji wielomianowej dla zestawu danych Głębia Challengera (węzły Czebyszewa)]
)
#pagebreak()
== Interpolacja funkcjami sklejanymi trzeciego stopnia

#figure(
  image("plots/chelm/cubic_spline.png"),
  caption: [Wpływ liczby punktów pomiarowych na dokładność interpolacji funkcjami sklejanymi dla zestawu danych Chelm]
)

#figure(
  image("plots/GlebiaChallengera/cubic_spline.png"),
  caption: [Wpływ liczby punktów pomiarowych na dokładność interpolacji funkcjami sklejanymi dla zestawu danych Głębia Challengera]
)
#pagebreak()
= Analiza dodatkowa interpolacji

W ramach analizy dodatkowej, przeprowadzono testy dla zestawu danych `100.csv`, ponieważ jest to zestaw danych o niezwykłym charakterze trasy, co pozwala na lepsze zbadanie zachowania obu metod interpolacji.
Przy interpolacji wielomianowej Lagrange'a widać wyraźny efekt Rungego jak w przypadku stosowania węzłów równoodległych (przy 23 punktach), tak i w przypadku stosowania węzłów Czebyszewa (przy 103 punktach). W przypadku interpolacji funkcjami sklejanymi trzeciego stopnia, efekt Rungego nie występuje. W każdej metodzie interpolacja jest niedokładna w niektórych przedziałach trasy z powodu skomplikowanego profilu terenu, jednak interpolacja funkcjami sklejanymi trzeciego stopnia daje najlepsze rezultaty.


#figure(
  image("plots/100/linspace_nodes.png"),
  caption: [Wpływ liczby punktów pomiarowych na dokładność interpolacji wielomianowej dla zestawu danych 100 (węzły równoodległe)]
)

#figure(
  image("plots/100/chebyshev_nodes.png"),
  caption: [Wpływ liczby punktów pomiarowych na dokładność interpolacji wielomianowej dla zestawu danych 100 (węzły Czebyszewa)]
)

#figure(
  image("plots/100/cubic_spline.png"),
  caption: [Wpływ liczby punktów pomiarowych na dokładność interpolacji funkcjami sklejanymi dla zestawu danych 100]
)

= Podsumowanie

Najmniej dokładna okazała się interpolacja wielomianowa Lagrange'a z równoodległymi węzłami. Efekt Rungego jest bardzo widoczny już przy stosunkowo niewielkiej liczbie punktów pomiarowych, co ogranicza skuteczność tej metody. Jest ona użyteczna jedynie w przypadku najprostszych funkcji terenu z małą liczbą punktów pomiarowych.

Interpolacja wielomianowa Lagrange'a z węzłami Czebyszewa przynosi zadowalające wyniki.
Pomaga uniknąć efektu Rungego do pewnej liczby węzłów interpolacji, jednak nie jest w stanie zapewnić takiej dokładności jak funkcje sklejane trzeciego stopnia.

Najlepsze wyniki daje interpolacja funkcjami sklejanymi trzeciego stopnia, zarówno pod względem dokładności, jak i odporności na efekt Rungego. Ta metoda zapewnia najwyższą jakość interpolacji w każdych warunkach. Chociaż może uprościć profil terenu w przypadku bardziej skomplikowanych funkcji, jest to zazwyczaj akceptowalne.
