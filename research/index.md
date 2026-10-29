---
layout: default
title: Research
---

# Research

<p class="lead">I use algebra, geometry and combinatorics to understand statistical models: what they look like, how to represent them, and when their estimators have closed-form solutions.</p>

<div class="theme" markdown="1">
<div markdown="1">
### Staged trees and context-specific models

Staged trees represent statistical models in which some independences hold only in certain contexts. I study their defining equations, when two trees describe the same model, and how to learn them from observational and interventional data.

<p class="meta"><b>With</b> {% include people.html names="Liam Solus, Christiane Görgen, Lamprini Ananiadi" %}<br>
<b>Key papers</b> <a href="https://arxiv.org/abs/2101.09271">Representation of context-specific causal models</a> · <a href="https://www.sciencedirect.com/science/article/abs/pii/S0747717119300379">Equations defining probability tree models</a> · <a href="https://arxiv.org/abs/1910.02721">Gröbner bases for staged trees</a></p>
</div>
<figure class="fig">
{% include staged-tree.svg %}
<figcaption>Nodes of the same colour share a stage.</figcaption>
</figure>
</div>

<div class="theme" markdown="1">
### Graphical and causal models

Graphical models encode dependence between random variables with a graph. I study their algebraic structure, from discrete models with interventions to Gaussian polytrees, and recently to quantum many-body systems.

<p class="meta"><b>With</b> {% include people.html names="Liam Solus, Yulia Alexandr, Julian Vill, Daniele Tramontano, Mathias Drton, Dmitrii Pavlov, Maximilian Wiesmann" %}<br>
<b>Key papers</b> <a href="https://arxiv.org/abs/2012.03593">Algebraic geometry of discrete interventional models</a> · <a href="https://ieeexplore.ieee.org/document/10299801">Learning linear Gaussian polytree models with interventions</a> · <a href="https://arxiv.org/abs/2308.11538">Algebraic geometry of quantum graphical models</a></p>
</div>

<div class="theme" markdown="1">
### Maximum likelihood and algebraic optimization

Which discrete statistical models have a maximum likelihood estimator given by a rational formula? This question connects statistics with toric geometry and polytopes. I also study related optimization problems, such as finding the closest point to a model in polyhedral norms.

<p class="meta"><b>With</b> {% include people.html names="Orlando Marigliano, Bernd Sturmfels, Isobel Davies, Irem Portakal, Miruna-Ştefana Sorea, Nidhi Kaihnsa, Julia Lindberg, Angélica Torres, Madeleine Weinstein" %}<br>
<b>Key papers</b> <a href="https://arxiv.org/abs/1903.06110">Discrete statistical models with rational MLE</a> · <a href="https://arxiv.org/abs/2109.08151">Polytopes with rational linear precision</a> · <a href="https://arxiv.org/abs/2604.19479">Algebraic distance optimization in polyhedral norms</a></p>
</div>

<div class="theme" markdown="1">
### Geometric modeling and implicitization

Curves and surfaces in computer-aided design are given by parametrizations. Finding their implicit equations uses syzygies and free resolutions from commutative algebra.

<p class="meta"><b>With</b> {% include people.html names="Hal Schenck, Alexandra Seceleanu, Benjamin Hollering, Maximilian Wiesmann" %}<br>
<b>Key papers</b> <a href="https://arxiv.org/abs/1402.6751">Tensor product surfaces and linear syzygies</a> · <a href="https://arxiv.org/abs/1908.02086">Implicitization via virtual projective resolutions</a> · <a href="https://arxiv.org/abs/2303.08754">Toric fiber products in geometric modeling</a></p>
</div>

## Software

<div class="cards">
  <div class="card"><h3><a href="https://github.com/emduart2/staged-tree">staged-tree</a></h3><p>A formalization of staged trees in the Lean 4 proof assistant.</p></div>
  <div class="card"><h3><a href="https://github.com/emduart2/DiscreteStatisticalModelsWithRationalMLE">Rational MLE</a></h3><p>Code accompanying <em>Discrete statistical models with rational MLE</em>.</p></div>
</div>
