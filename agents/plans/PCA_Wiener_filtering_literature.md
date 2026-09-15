I'm studying Wiener filtering in the context of high contrast imaging post processing.  We have PCA of our data (sets of images).  There is a deep connection between PCA/eigendecomposition and PSD estimation, so it seems to me that since I can easily produce the projection of both the noise and the signal onto the PCA basis, I should be able to use that to construct the Wiener filter rather than using the Fourier-basis PSD as it is classically described.  Is there any literature on this?

Yes. What you are describing is a well-established generalization of Wiener filtering, usually formulated as Wiener filtering in the Karhunen–Loève (KL) domain or as covariance-domain Wiener filtering. PCA gives you the empirical KL basis. The Fourier-domain version is really a special case that becomes convenient for stationary processes.

The key point is that the Wiener filter fundamentally depends on covariances, not Fourier PSDs.

Suppose your data vectorized into images is

$$ \mathbf y = \mathbf s + \mathbf n, $$

with signal and noise covariances

$$ C_s = E[\mathbf s\mathbf s^\mathsf{T}],\qquad C_n = E[\mathbf n\mathbf n^\mathsf{T}]. $$

The minimum-MSE linear estimator is

$$ \hat{\mathbf s} = W\mathbf y, \qquad W=C_s(C_s+C_n)^{-1}. $$

This is the generalized Wiener filter. The familiar Fourier expression

$$ W(f)=\frac{P_s(f)}{P_s(f)+P_n(f)} $$

comes from the fact that, for stationary processes, the Fourier basis diagonalizes (asymptotically, in the finite discrete case) the covariance operators. So the PSDs are essentially just the covariance eigenvalues expressed in that particular basis. Pratt's classic 1972 paper explicitly formulates generalized Wiener filtering using arbitrary unitary transforms, including Fourier, Hadamard, and Karhunen–Loève transforms.

Your PCA construction

Let the PCA/KL eigenvectors be the columns of

$$ U=[\mathbf u_1,\ldots,\mathbf u_N]. $$

Transform the signal and noise:

$$ \mathbf a_s=U^\mathsf{T}\mathbf s,\qquad \mathbf a_n=U^\mathsf{T}\mathbf n. $$

Then

$$ C_s' = U^\mathsf{T}C_sU,\qquad C_n' = U^\mathsf{T}C_nU. $$

The Wiener estimator in this basis is simply

$$ \hat{\mathbf a}_s = C_s'(C_s'+C_n')^{-1}\mathbf a_y, $$

followed by

$$ \hat{\mathbf s}=U\hat{\mathbf a}_s. $$

So you do not actually need a PSD at all.

There is one important qualification to your idea. If the PCA basis simultaneously diagonalizes the signal and noise covariances, then this becomes exactly the scalar operation you have in mind:

$$ \hat a_{s,k} = \frac{\sigma_{s,k}^2} {\sigma_{s,k}^2+\sigma_{n,k}^2} a_{y,k}. $$

That is literally the ordinary Wiener gain with

$$ P_s(f)\longrightarrow \sigma_{s,k}^2, \qquad P_n(f)\longrightarrow \sigma_{n,k}^2. $$

Chen, Benesty & Huang (2009) develop essentially exactly this formulation in the KL expansion domain: transform into KL components, estimate the signal/noise power in each component, Wiener-filter the components, and transform back. They explicitly contrast this with the Fourier-domain formulation and note that the KLT exactly diagonalizes the relevant correlation matrix whereas the Fourier transform generally only approximately does so.

That paper is probably the single best starting reference for the mathematical idea you have in mind:

Chen, Benesty & Huang 2009 — Study of the Noise-Reduction Problem in the Karhunen–Loève Expansion Domain

There is also an earlier/general treatment by Hua & Liu, whose generalized KL transform explicitly connects rank-reduced KL estimation to Wiener filtering.

The subtle issue for your high-contrast application

I think this gets particularly interesting in your case because you presumably construct PCA primarily from the speckle/noise ensemble.

If

$$ C_n = U\Lambda_nU^\mathsf{T}, $$

then

$$ C_n'=\Lambda_n $$

is diagonal by construction. But generally

$$ C_s'=U^\mathsf{T}C_sU $$

will not be diagonal.

If you merely calculate the planet/signal projection onto each PCA mode and assign a signal power

$$ P_{s,k}=E[|a_{s,k}|^2], $$

then

$$ g_k = \frac{P_{s,k}} {P_{s,k}+\lambda_{n,k}} $$

is a diagonal approximation to the full Wiener filter.

That may nevertheless be exactly what you want computationally. It says:

retain each empirical speckle eigenmode according to the amount of astrophysical signal relative to speckle variance carried by that mode.

That is substantially more nuanced than ordinary PCA truncation. PCA/KLIP effectively makes something resembling a hard decision about modes, whereas this gives you a continuous shrinkage:

$$ 0\le g_k\le1. $$

There's a nice conceptual progression:

$$ \text{KLIP/PCA truncation} \quad\longrightarrow\quad \text{KL-mode shrinkage} \quad\longrightarrow\quad \text{full covariance Wiener estimator}. $$
Connection to existing high-contrast work

I did not immediately find a high-contrast-imaging paper that says, in precisely these terms, "take the KLIP basis and apply Wiener gains derived from planet and speckle projections." That may be why the idea feels simultaneously obvious and unfamiliar.

But several pieces of the HCI literature sit very close to it.

Soummer, Pueyo & Larkin's original KLIP formulation is, of course, explicitly a KL decomposition of the empirical speckle covariance. Pueyo's KLIP forward-modeling work then calculates how an astrophysical signal projects through and perturbs that KL representation. Those give you essentially the ingredients on the PCA side.

More interestingly, PACO moves toward the other side of the same statistical picture. Flasseur et al. explicitly estimate local speckle covariance matrices and use those covariance models for optimal statistical detection rather than simply subtracting PCA modes. PACO-ASDI extends this to spatial/spectral/temporal covariance. In other words, PACO is closely related to what happens when you stop throwing away the off-diagonal covariance information.

There is also a striking analogue outside astronomy: Bhamre, Zhang & Singer call their cryo-EM method Covariance Wiener Filtering (CWF). They estimate the covariance of an ensemble of images and insert it directly into the Wiener estimator rather than relying on a conventional Fourier PSD. Conceptually that is extremely close to what you're proposing.

One distinction I'd make for planet detection

If by "signal" you mean a known planet PSF/template at a candidate location, then I would be slightly careful about calling its PCA projections a signal PSD. A deterministic planet PSF doesn't naturally have a covariance unless you define an ensemble—locations, parallactic angles, spectra, flux distributions, etc.

For detection of a known template \(\mathbf p\) in correlated speckle noise, the natural statistic is instead the covariance-weighted matched filter,

$$ T = \frac{\mathbf p^\mathsf{T}C_n^{-1}\mathbf y} {\sqrt{\mathbf p^\mathsf{T}C_n^{-1}\mathbf p}}. $$

In your PCA basis,

$$ C_n=U\Lambda_nU^\mathsf{T}, $$

so this becomes

$$ T= \frac{\displaystyle \sum_k \frac{p_k y_k}{\lambda_k}} {\displaystyle \sqrt{\sum_k\frac{p_k^2}{\lambda_k}}}, $$

where

$$ p_k=\mathbf u_k^\mathsf{T}\mathbf p,\qquad y_k=\mathbf u_k^\mathsf{T}\mathbf y. $$

That is an especially attractive formulation for HCI because you already have all of those quantities from the PCA. And this is closely connected to the covariance-aware detection philosophy in PACO.

If your goal instead is reconstruction/estimation of the astrophysical image, rather than detection of a known planet template, then the KL-domain Wiener filter

$$ g_k= \frac{\sigma_{s,k}^2} {\sigma_{s,k}^2+\lambda_{n,k}} $$

is exactly the natural object.

So I think your intuition is correct, with the important distinction:

$$ \boxed{\text{Fourier PSD Wiener filtering} \;\;\text{is a special case of covariance/KL-domain Wiener filtering}.} $$

And for high-contrast imaging, the interesting question isn't whether this is mathematically valid—it definitely is—but rather what signal covariance you want to define, and whether the diagonal approximation in the speckle PCA basis is adequate.

The Chen et al. 2009 paper is where I'd start mathematically; then I'd put it next to KLIP/KLIP-FM and PACO. That triangle gets very close to the exact method you're contemplating.

Since you're collecting literature around a fairly specific methodological connection, Zotero could also be useful for keeping these papers and citations together.

Wed, Aug 26 at 11:51 AM
in the covar-weighted matchef filter in the PCA basis, the u\_k are eigenvectors and lambda\_k the eigenvalues?

Yes, exactly.

If the noise/speckle covariance matrix has eigendecomposition

$$ C_n = U\Lambda U^\mathsf{T}, $$

then

$$ U = [\mathbf u_1,\mathbf u_2,\ldots] $$

contains the PCA eigenvectors, and

$$ \Lambda = \operatorname{diag}(\lambda_1,\lambda_2,\ldots) $$

contains the corresponding eigenvalues. Thus

$$ C_n^{-1} = U\Lambda^{-1}U^\mathsf{T}. $$

For a planet template \(\mathbf p\) and observed image \(\mathbf y\), define their PCA coefficients

$$ p_k = \mathbf u_k^\mathsf{T}\mathbf p, \qquad y_k = \mathbf u_k^\mathsf{T}\mathbf y. $$

Then

$$ \mathbf p^\mathsf{T}C_n^{-1}\mathbf y = \sum_k \frac{p_k y_k}{\lambda_k}, $$

and the normalized matched-filter statistic is

$$ T = \frac{ \displaystyle\sum_k \frac{p_k y_k}{\lambda_k} }{ \displaystyle\sqrt{\sum_k\frac{p_k^2}{\lambda_k}} }. $$

So each PCA mode is weighted by inverse noise variance \(1/\lambda_k\).

This makes the connection particularly clean: the PCA eigenvalue \(\lambda_k\) is the noise power in mode \(k\), playing essentially the same role that the noise PSD \(P_n(f)\) plays for Fourier modes.

One caveat: this assumes the PCA basis was obtained from the noise covariance \(C_n\), rather than from signal+noise data in which the astrophysical signal significantly contributes to the covariance.