Landmarking Protocol 3d3 + Hypothesis Viz
================
Robert Z. Selden, Jr.
05 September, 2020

Landmarking protocol 3d3 (LM3d3) represents a substantive advancement
from
[LM3d1](https://github.com/aksel-blaise/gahaganmorph2/blob/master/analysis/landmarking-protocol.md)
(Selden Jr., Dockall, and Dubied 2020) and was modified from the
[LM3d2](https://aksel-blaise.github.io/gahaganmorph3/landmarking-protocol.html)
protocol developed for the analysis of Gahagan bifaces. The principal
difference between
[LM3d2](https://aksel-blaise.github.io/gahaganmorph3/landmarking-protocol.html)
and LM3d3 is that the projectiles analysed using
[LM3d2](https://aksel-blaise.github.io/gahaganmorph3/landmarking-protocol.html)
are lanceolate bifaces, and that landmarking protocol was not designed
to capture the variation in basal morphology exhibited by stemmed
(Dalton) projectile points. Like
[LM3d1](https://github.com/aksel-blaise/gahaganmorph2/blob/master/analysis/landmarking-protocol.md)
and
[LM3d2](https://aksel-blaise.github.io/gahaganmorph3/landmarking-protocol.html),
LM3d3 uses the topology of the 3D mesh that articulates with the
prehistoric design of each projectile point to construct a suite of
`reference geometry` used to apply semilandmarks in a replicable manner.
The result is a landmarking protocol that provides for the improved
characterisation of whole-object morphology, which can be subset to
analyse variability associated with specific *plan*, *profile*, and
*cross-section* components. It also provides those data points needed to
investigate questions of *directional asymmetry*, differences in
front/back morphology, morphological integration between blade and basal
morphology, and the morphology of broken, fractured, or otherwise
incomplete specimens that permeate the archaeological record.

``` r
knitr::include_graphics('images/landmarks.png')
```

<div class="figure">

<img src="images/landmarks.png" alt="_**Figure 1. Coordinates of landmarks populated using the [LM3d1](https://github.com/aksel-blaise/gahaganmorph2/blob/master/analysis/landmarking-protocol.md) protocol serve as the basis for LM3d3.**_" width="100%" />

<p class="caption">

***Figure 1. Coordinates of landmarks populated using the
[LM3d1](https://github.com/aksel-blaise/gahaganmorph2/blob/master/analysis/landmarking-protocol.md)
protocol serve as the basis for LM3d3.***

</p>

</div>

Like its’ predecessors, LM3d3 was initially designed using the
[`digit3DLand`](https://github.com/morphOptics/digit3DLand) package in
R. When the draft protocol was completed, design of the landmarking
protocol shifted to [Geomagic Design
X](https://www.3dsystems.com/software/geomagic-design-x) *(Build Version
2020.0.1 \[Build Number: 30\])*, where the workflow was modified to
include those elements of `reference geometry` that articulate with the
prehistoric design attributes of each projectile point.

The goal of this effort was to increase both the precision and rigour of
the study by including the Z-dimension to capture those shape
characteristics associated with axial twisting, introduced by knappers
through the practice of beveling
([LM3d1](https://github.com/aksel-blaise/gahaganmorph2/blob/master/analysis/landmarking-protocol.md))
(Selden Jr., Dockall, and Dubied 2020). The addition of cross-sections
was needed to better characterise whole-object morphology, providing for
the possibility of subsampling the semilandmarks to explore the
contribution of specific cross-sections or profiles
([LM3d2](https://aksel-blaise.github.io/gahaganmorph3/landmarking-protocol.html)).
LM3d3 includes an additional cross section at the blade/base transition,
allowing for tests of morphological integration. While true that some
landmarking protocols can be—–and often are—–recycled as new specimens
are added, this particular research programme endeavours to achieve
ever-greater accuracy and precision in each analytical iteration.

## Generating the peripheral (plan view) spline

This effort begins with a spline extracted from the surface geometry of
the mesh using the `extract contour curves` command. In
reverse-engineering, `extract contour curves` is regularly employed as
the first step in building a `patch network` used to create a surface.
The extracted feature curve is rendered as a spline, and follows the
highest curvature contours around the periphery of the lateral and basal
edges, following the highly variable sinuous edge morphology around the
entirety of the projectile. The remainder of the landmarking protocol is
based upon this spline, which was subsequently split at six
mathematically-defined locations.

``` r
knitr::include_graphics('images/extractspline.png')
```

<div class="figure">

<img src="images/extractspline.png" alt="_**Figure 2. Spline extracted along the highest contours of the Dalton point.**_" width="100%" />

<p class="caption">

***Figure 2. Spline extracted along the highest contours of the Dalton
point.***

</p>

</div>

## Splitting the spline

*`Reference geometries` are used in the assistance of creating other
features. These include basic geometric entities, such as `planes`,
`vectors`, `coordinates`, `points`, and `polygons`. A `reference point`
is a virtual point and is used to mark a specific position on a model or
in 3D space. A `reference plane` is a virtual plane that has a normal
direction and an infinite size. A `reference plane` is not a surface
body, and is used to create other features.*

The characteristic points and tangents developed for this landmarking
protocol were inspired by the work of Birkhoff (1933), which has been
gainfully employed within the context of both ceramic (Selden Jr. 2018a,
2018b, 2019, 2021) and lithic analyses (Selden Jr., Dockall, and Shafer
2018; Selden Jr., Dockall, and Dubied 2020). The first landmark (LM1) is
placed at the horizontal tangent on the tip of each Dalton point. The
second through fifth splits (LMs 02 - 05) occur at points of highest
curvature, where LM 02 is always placed on the right side of the
projectile following the application of the `reference vectors`. To
place the final landmark (LM 06), a linear measurement was used to
project a `reference point` equidistant between LM 02 and LM 03. The
location of that point was leveraged in placing the `reference plane`
used to cut the spline at the location of LM 06.

## Spline split at location of LM 01

The `horizontal tangent` is calculated by drawing a horizontal line
above the tip of the biface using the tangent as a `common constraint`,
and the horizontal as the `independent constraint`. To split the 3D
spline at the location of the horizontal tangent, a `reference point`
was inserted at the location of the `tangent` in the sketch (light blue
point; below, left), followed by a `reference plane` (in white; below,
left and right) using the `pick point and normal axis` function where
the `reference point` (h-tangent) was used as the `pick point`, and the
`Right plane` as the `normal axis` (below, left). The spline was then
cut at the location where the `reference plane` intersected with the
spline (below image, right).

``` r
knitr::include_graphics('images/lm1.png')
```

<div class="figure">

<img src="images/lm1.png" alt="_**Figure 3. Identify horizontal tangent, insert `reference point` and `reference plane` (left). Use `reference plane` to cut spline at the location of the horizontal tangent (right).**_" width="100%" />

<p class="caption">

***Figure 3. Identify horizontal tangent, insert `reference point` and
`reference plane` (left). Use `reference plane` to cut spline at the
location of the horizontal tangent (right).***

</p>

</div>

## Spline split at locations of LM 02 and LM 03

The point of highest curvature on either side of the basal edge was
calculated using the `curvature` function in the Accuracy Analyser. This
function displays the curvature flow as a continuous colour plot across
the area of the curve. In this instance, *curvature* is defined as the
amount by which a geometric shape deviates from being flat or straight
in the case of a line. Curvature is displayed in different colours
according to the local radius, and calculated in only one direction (U
or V) along the curve. Using this tool, the two points of highest
curvature were located between the basal and lateral edges on either
side of each projectile where the local radius measure was largest. The
orientation of each biface was dictated by the *auto3dgm* output in
[LM3d1](https://github.com/aksel-blaise/gahaganmorph2/blob/master/analysis/landmarking-protocol.md)
and
[LM3d2](https://aksel-blaise.github.io/gahaganmorph3/landmarking-protocol.html);
however, LM3d3 enlists a novel method to determine which side of the
projectile is associated with LM 02 and LM03 using `reference vectors`.

``` r
knitr::include_graphics('images/splinesplit1.png')
```

<div class="figure">

<img src="images/splinesplit1.png" alt="_**Figure 4. Identify points of hightest curvature (light blue) at left/right intersection of lateral and basal edges.**_" width="100%" />

<p class="caption">

***Figure 4. Identify points of hightest curvature (light blue) at
left/right intersection of lateral and basal edges.***

</p>

</div>

## Spline split at locations of LM 04 and LM 05

The point of highest curvature at the intersection of the blade and base
was also calculated using the `curvature` function in the Accuracy
Analyser. Using this tool, the two points of highest curvature were
located between the blade and base on either side of each projectile
where the local radius measure was largest. The orientation of each
projectile was dictated by `reference vectors`, and the landmarking
protocol follows the mesh orientation in that figure, where LM 04 was
always placed on the right side of the basal edge, and LM 05 on the
left.

``` r
knitr::include_graphics('images/splinesplit2.png')
```

<div class="figure">

<img src="images/splinesplit2.png" alt="_**Figure 5. Identify points of hightest curvature (light blue) at left/right intersection of blade and base.**_" width="100%" />

<p class="caption">

***Figure 5. Identify points of hightest curvature (light blue) at
left/right intersection of blade and base.***

</p>

</div>

## Spline split at location of LM 06

One additional landmark (LM 06) was placed at the centre of the base.
The location of this landmark was identified by calculating the linear
distance between LM 02 and LM 03, and projecting a `reference point`
(ctrl-div; below) equidistant between the two. A `reference plane` was
added using the ctrl-div as the pick point, and the `Right plane` as the
`normal axis`. The spline was then split at the intersection of the
`reference plane` and the basal spline.

``` r
knitr::include_graphics('images/lm6.png')
```

<div class="figure">

<img src="images/lm6.png" alt="_**Figure 6. Calculate linear distance between LM 02 and LM 03, insert `reference plane` coplanar to Right plane equidistant between LM 02 and LM 03, and use the `reference plane` to cut the spline.**_" width="100%" />

<p class="caption">

***Figure 6. Calculate linear distance between LM 02 and LM 03, insert
`reference plane` coplanar to Right plane equidistant between LM 02 and
LM 03, and use the `reference plane` to cut the spline.***

</p>

</div>

## Peripheral (plan view) spline

Through the preceding protocol, the initial spline was split into six
discrete splines. These splines articulate with components of projectile
morphology that can be compartmentalised in the analyses. The primary
analytical gain achieved through this exercise is the requisite
foundation needed to carry out replicable analyses of Dalton point
morphology in three dimensions, further increasing the precision of the
geometric morphometric analysis.

``` r
knitr::include_graphics('images/splinesplit-frbl.png')
```

<div class="figure">

<img src="images/splinesplit-frbl.png" alt="_**Figure 7. Result of spline splits include six discrete splines, each articulating with a region of analytical interest. The coordinates of each spline split are known, and used to place the landmarks.**_" width="100%" />

<p class="caption">

***Figure 7. Result of spline splits include six discrete splines, each
articulating with a region of analytical interest. The coordinates of
each spline split are known, and used to place the landmarks.***

</p>

</div>

## Reference vectors and ref.pt.0

The fundamental components of `reference geometry` used to create LM3d3
consist of three `reference vectors`, and a single `reference point`
(ref.pt.0), placed equidistant between LM 04 and LM 05. The three
`reference vectors` were placed between LM 01 and ref.pt.0 (Vector 1),
ref.pt.0 and LM 06 (Vector 2), and LMs 04 and 05 (Vector 3). These three
`reference vectors` serve as the foundation for the suite of `reference
geometry` used to place the semilandmarks.

``` r
knitr::include_graphics('images/lm3d3.vectors.png')
```

<div class="figure">

<img src="images/lm3d3.vectors.png" alt="_**Figure 8. `Reference vectors` placed between LMs 01 and ref.pt.0 (left), ref.pt.0 and LM 06 (center), and between LMs 04 and 05 (right).**_" width="100%" />

<p class="caption">

***Figure 8. `Reference vectors` placed between LMs 01 and ref.pt.0
(left), ref.pt.0 and LM 06 (center), and between LMs 04 and 05
(right).***

</p>

</div>

The measure of the angle between Vector 1 (blade) and Vector 2 (base)
may have additional utility in lithic studies as an orthogonal metric
associated with knapper skill, where greater skill is represented by an
arbitrary range of angles nearest—and lesser, furthest away from—180
degrees. A second similar measure could be collected between Vectors 1
and 2 in comparison with Vector 3. Collection of these metrics from a 3D
mesh in computer aided design (CAD) software adds an increased element
of precision in comparison with a goniometer, and serves as an example
of the added analytical value that can be extracted from this novel
landmarking protocol.

Prior to the addition of the `reference vectors`, the location of LMs 02
through 05 are considered arbitrary. Previous iterations of this
landmarking protocol have relied upon `auto3dgm` to provide principal
alignments that dictate which LMs are placed on the left/right side of
the biface or projectile. In this protocol, the side of the projectile
with the lowest orthogonal measure between Vector 1 and Vector 2 will be
on the left, meaning that from the investigator’s view, the projectiles
will curve, bend, or lean slightly—or in some cases more dramatically—to
the left from base to tip.

## Reference planes and points

Five `reference planes` provide the framework needed to populate the
semilandmarks. Admittedly, the logic associated with placement may seem
curious at this point; however, the utility of these `reference planes`
will become clear in subsequent sections.

### Placement of ref.pl.1

The first `reference plane` (ref.pl.1) was placed between LM 01 and
ref.pt.0, bisecting the blade of the projectile along the mid-line. The
method of placement enlists a second `reference point` (ref.pt.1),
inserted along the first `reference vector`. It is located at a position
equidistant between LM 01 and ref.pt.0, but the coordinates of ref.pt.1
were altered to relocate it 15 mm from the vector in the direction of
the Z-axis. The `pick point and coplanar` function was used to place
ref.pl.1 coplanar to the first `reference vector`, and in the direction
of ref.pt.1. Following placement of ref.pl.1, ref.pt.1 was deleted.

``` r
knitr::include_graphics('images/lm3d3.ref.pl.1.png')
```

<div class="figure">

<img src="images/lm3d3.ref.pl.1.png" alt="_**Figure 9. Placement of ref.pl.1, and temporary location of ref.pt.1 15mm from Vector 1, and equidistant between LM 01 and ref.pt.0 on the blade of the projectile.**_" width="100%" />

<p class="caption">

***Figure 9. Placement of ref.pl.1, and temporary location of ref.pt.1
15mm from Vector 1, and equidistant between LM 01 and ref.pt.0 on the
blade of the projectile.***

</p>

</div>

### Placement of ref.pl.2

The second `reference plane` (ref.pl.2) was placed between ref.pt.0 and
LM 06, bisecting the base of the projectile along the mid-line. The
method of placement for ref.pl.2 follows the same protocol described in
the application of ref.pl.1, and the `reference point` (ref.pt.2) was
deleted following the placement of ref.pl.2.

``` r
knitr::include_graphics('images/lm3d3.ref.pl.2.png')
```

<div class="figure">

<img src="images/lm3d3.ref.pl.2.png" alt="_**Figure 10. Placement of ref.pl.2, and temporary location of ref.pt.2 15 mm from Vector 2, and equidistant between ref.pt.0 and LM 06 on the base of the projectile.**_" width="100%" />

<p class="caption">

***Figure 10. Placement of ref.pl.2, and temporary location of ref.pt.2
15 mm from Vector 2, and equidistant between ref.pt.0 and LM 06 on the
base of the projectile.***

</p>

</div>

### Placement of ref.pl.3

The third `reference plane` (ref.pl.3) was placed between LMs 04 and 05,
and bisects the projectile at the blade/base intersection. The method of
placement for ref.pl.3 follows the same protocol described in the
application of ref.pl.1, and the `reference point` (ref.pt.3) was
deleted following the placement of ref.pl.3.

``` r
knitr::include_graphics('images/lm3d3.ref.pl.3.png')
```

<div class="figure">

<img src="images/lm3d3.ref.pl.3.png" alt="_**Figure 11. Placement of ref.pl.3, and temporary location of ref.pt.3 15 mm from Vector 3, and equidistant between LMs 04 and 05 at the intersection of the blade and base.**_" width="100%" />

<p class="caption">

***Figure 11. Placement of ref.pl.3, and temporary location of ref.pt.3
15 mm from Vector 3, and equidistant between LMs 04 and 05 at the
intersection of the blade and base.***

</p>

</div>

### Placement of ref.pl.4 and ref.pl.5

The fourth (ref.pl.4) and fifth (ref.pl.5) `reference planes` were
placed using the `pick point and normal` function at the intersections
of the first `reference vector` and LMs 01 and 06.

``` r
knitr::include_graphics('images/lm3d3.ref.pl.4-5.png')
```

<div class="figure">

<img src="images/lm3d3.ref.pl.4-5.png" alt="_**Figure 12. Placement of ref.pl.4 (top) and ref.pl.5 (bottom).**_" width="100%" />

<p class="caption">

***Figure 12. Placement of ref.pl.4 (top) and ref.pl.5 (bottom).***

</p>

</div>

## Sectioning the mesh

The `reference geometry` described above was enlisted in the following
three-step method developed to produce one cross-section at the
blade/base intersection, four cross-sections between the blade/base
intersection and LM 01, and one cross-section between the blade/base
intersection and LM 06.

### Sectioning the blade/base intersection

To section the blade/base intersection, a single section was inserted
using ref.pl.2, resulting in a single cross-section that bisects the
projectile between LMs 04 and 05.

``` r
knitr::include_graphics('images/lm3d3.section1.png')
```

<div class="figure">

<img src="images/lm3d3.section1.png" alt="_**Figure 13. Placement of the first section, bisecting the mesh along ref.pl.2.**_" width="100%" />

<p class="caption">

***Figure 13. Placement of the first section, bisecting the mesh along
ref.pl.2.***

</p>

</div>

### Sectioning the blade and base

Six equidistant sections were placed between LM 01 and ref.pt.0, and the
two sections at the locations of LM 01 and ref.pt.0 were deleted. Three
equidistant sections were placed between LM 06 and ref.pt.0. The
sections intersecting with ref.pt.0 and LM 06 were deleted. Subsequent
to placing the sections, ref.pt.0 was itself deleted.

``` r
knitr::include_graphics('images/lm3d3.all.sections.png')
```

<div class="figure">

<img src="images/lm3d3.all.sections.png" alt="_**Figure 14. Placement of the two equidistant sections between LM 06 and ref.pt.3.**_" width="100%" />

<p class="caption">

***Figure 14. Placement of the two equidistant sections between LM 06
and ref.pt.3.***

</p>

</div>

## Splitting the sections

The `curvature` function was employed to split each curves at the
locations of highest curvature along the lateral edge. This function was
detailed above, and in the application of LMs 02, 03, 04, and 05 in
[LM3d1](https://github.com/aksel-blaise/gahaganmorph2/blob/master/analysis/landmarking-protocol.md).

A `reference plane` (ref.pl.1) was then used to cut each of the four
curves along the mid-line of the blade where it intersects with the
curves. A second `reference plane` (ref.pl.2) was used to cut the single
basal curve and the curve between LMs 04 and 05. Since ref.pt.0 was used
to generate ref.pl.1 and ref.pl.2, either of the `reference planes`
could be used to cut the curve between LMs 04 and 05.

``` r
knitr::include_graphics('images/lm3d3.split.sections.png')
```

<div class="figure">

<img src="images/lm3d3.split.sections.png" alt="_**Figure 15. Each section was split at the points of highest curvature along the lateral edges, then along the mid-line at the intersection of the curve and ref.pl.1 (for the blade), and ref.pl.2 (for the base).**_" width="100%" />

<p class="caption">

***Figure 15. Each section was split at the points of highest curvature
along the lateral edges, then along the mid-line at the intersection of
the curve and ref.pl.1 (for the blade), and ref.pl.2 (for the base).***

</p>

</div>

## LM3d3: Configuration 1

LM3d3: Configuration 1 (LM3d3:c1) was used to assess the first three
hypotheses. Semilandmarks 07 - 18 were first applied around the lateral
edges, and sLMs 12 and 13 between LMs 02, 06, and 03 uses the curve
constructed in
[LM3d1](https://github.com/aksel-blaise/gahaganmorph2/blob/master/analysis/landmarking-protocol.md).
Semilandmarks 19 - 30 articulate with the mid-line between LMs 01 and
06. Additional sLMs can be placed on each section equidistant between
those sLMs defined by the splits described above to better characterise
each cross section, if needed.

The result is a landmark configuration that can be subset in numerous
ways (plan, profile, cross-section, front/back, left/right, blade/base,
etc.), and was designed to achieve maximum utility for analysts of
lithic morphology.

``` r
knitr::include_graphics('images/lm3d3.semi.png')
```

<div class="figure">

<img src="images/lm3d3.semi.png" alt="_**Figure 16. Landmarks (blue), semilandmarks (white), curves (orange), and splits (blue) used for LM3d3.**_" width="100%" />

<p class="caption">

***Figure 16. Landmarks (blue), semilandmarks (white), curves (orange),
and splits (blue) used for LM3d3.***

</p>

</div>

Based upon knowledge garnered from running `LaSEC` (Watanabe 2018) on
[LM3d1](https://github.com/aksel-blaise/gahaganmorph2/blob/master/analysis/landmarking-protocol.md),
this landmarking protocol would likely be oversampled if it included
additional landmarks on the cross-sections; however, it can be adapted
to include as many or as few landmarks and semilandmarks needed to
address the research question.

``` r
knitr::include_graphics('images/lm3d3.slm.png')
```

<div class="figure">

<img src="images/lm3d3.slm.png" alt="_**Figure 17. Landmarks (blue), semilandmarks (white), and `reference geometry` used in LM3d3.**_" width="100%" />

<p class="caption">

***Figure 17. Landmarks (blue), semilandmarks (white), and `reference
geometry` used in LM3d3.***

</p>

</div>

## Hypotheses 1 through 3

A suite of visual aids were created as a means of visualising the
hypotheses that LM3d3:c1 will be used to test. The visual representation
of the following hypotheses were produced for my own reference, and were
used to critically assess the utility of the semilandmarks employed in
the landmarking protocol. *All hypotheses tested for this study will
enlist the same sample of Dalton projectile points.*

### Hypothesis 1

Hypothesis 1 will test whether there is a difference in morphology for
Dalton points found in and out of the heartland.

``` r
knitr::include_graphics('images/dalton-vizhypothesis1.jpg')
```

<div class="figure">

<img src="images/dalton-vizhypothesis1.jpg" alt="_**Figure 18. Hypothesis 1 considers whether Dalton points discovered in (left) and out (right) of the Heartland differ in morphology.**_" width="100%" />

<p class="caption">

***Figure 18. Hypothesis 1 considers whether Dalton points discovered in
(left) and out (right) of the Heartland differ in morphology.***

</p>

</div>

### Hypothesis 2

Hypothesis 2 will test whether there is a difference in morphology for
Dalton points found in the heartland, the interior, and the northern
periphery.

``` r
knitr::include_graphics('images/dalton-vizhypothesis2.jpg')
```

<div class="figure">

<img src="images/dalton-vizhypothesis2.jpg" alt="_**Figure 19. Hypothesis 2 considers whether Dalton points discovered in the heartland (left), interior (center), and northern periphery (right) differ in morphology.**_" width="100%" />

<p class="caption">

***Figure 19. Hypothesis 2 considers whether Dalton points discovered in
the heartland (left), interior (center), and northern periphery (right)
differ in morphology.***

</p>

</div>

### Hypothesis 3

Hypothesis 3 will test whether there is a discernible difference in
morphology for Dalton points that are beveled, or not beveled.

``` r
knitr::include_graphics('images/dalton-vizhypothesis3.jpg')
```

<div class="figure">

<img src="images/dalton-vizhypothesis3.jpg" alt="_**Figure 20. Hypothesis 3 considers whether beveled Dalton points (left) differ in morphology from those that are not beveled (right).**_" width="100%" />

<p class="caption">

***Figure 20. Hypothesis 3 considers whether beveled Dalton points
(left) differ in morphology from those that are not beveled (right).***

</p>

</div>

## ***(in development)*** LM3d3: Configuration 2

LM3d3: Configuration 2 (LM3d3:c2) is an extension and subset of
LM3d3:c1, and is used in the assessment of Hypothesis 4.

## ***(in development)*** Sectioning the tip

Three additional equidistant sections were placed between LM 01 and the
first cross section on the blade, and were inserted using ref.pl.4 as
the base plane.

## ***(in development)*** Hypothesis 4

## Acknowledgments

I extend my gratitude to Christian S. Hoggard and David K. Thulman for
their thoughtful comments and constructive criticisms on the draft of
this landmarking protocol, which was originally developed for the study
of [Gahagan
bifaces](https://github.com/aksel-blaise/gahaganmorph2/blob/master/analysis/landmarking-protocol.md)
(Selden Jr., Dockall, and Dubied 2020), and is extended here for an
analysis of Dalton point morphology. This iteration of the landmarking
protocol was developed using the
[`digit3DLand`](https://github.com/morphOptics/digit3DLand) package in
R.

## References

<div id="refs" class="references">

<div id="ref-RN11786">

Birkhoff, George D. 1933. *Aesthetic Measure*. Cambridge: Harvard
University Press.

</div>

<div id="ref-RN11801">

Selden Jr., Robert Z. 2018a. “A Preliminary Study of Smithport Plain
Bottle Morphology in the Southern Caddo Area.” *Bulletin of the Texas
Archeological Society* 89: 63–89.

</div>

<div id="ref-RN11782">

———. 2018b. “Ceramic Morphological Organisation in the Southern Caddo
Area: Quiddity of Shape for Hickory Engraved Bottles.” *Journal of
Archaeological Science: Reports* 21: 884–96.
<https://doi.org/10.1016/j.jasrep.2018.08.045>.

</div>

<div id="ref-RN11716">

———. 2019. “Ceramic Morphological Organisation in the Southern Caddo
Area: The Clarence H. Webb Collections.” *Journal of Cultural Heritage*
35: 41–55.
<https://doi.org/https://doi.org/10.1016/j.culher.2018.07.002>.

</div>

<div id="ref-RN20697">

———. 2021. “Louisiana Limitrophe: An Iterative Morphological Exegesis of
Caddo Bottle and Biface Production.” In *Ancestral Caddo Ceramic
Traditions*, edited by Duncan P. McKinnon, Jeffrey S. Girard, and
Timothy K. Perttula, (in press). Baton Rouge: LSU Press.

</div>

<div id="ref-RN21001">

Selden Jr., Robert Z., John E. Dockall, and Morgane Dubied. 2020. “A
Quantitative Assessment of Intraspecific Morphological Variation in
Gahagan Bifaces from the Southern Caddo Area and Central Texas.”
*Southeastern Archaeology* 39 (2): 125–45.
<https://doi.org/10.1080/0734578x.2020.1744416>.

</div>

<div id="ref-RN11783">

Selden Jr., Robert Z., John E. Dockall, and Harry J. Shafer. 2018.
“Lithic Morphological Organisation: Gahagan Bifaces from the Southern
Caddo Area.” *Digital Applications in Archaeology and Cultural Heritage*
10: e00080. <https://doi.org/10.1016/j.daach.2018.e00080>.

</div>

<div id="ref-RN28913">

Watanabe, Akinobu. 2018. “How Many Landmarks Are Enough to Characterize
Shape and Size Variation?” Journal Article. *PLoS One* 13 (6): e0198341.
<https://doi.org/10.1371/journal.pone.0198341>.

</div>

</div>
