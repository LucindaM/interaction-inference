From: Sourya Shrestha <sshres14@jhu.edu>
To: Nicholas Reich <nick@schoolph.umass.edu>
Subject: Re: code
Date: Fri, 18 Jul 2014 18:43:48 +0000

Hi Nick,

I had meant to try to organize/comment this a little better, as well as make sure they run with newer pomp versions, but that might take a while.
But for now, I am just attaching the files that I had originally run to do one of the scenarios, may have been scenario I (no cross-serotype interaction).
I will try to get this compatible when I get a chance, but you may be able to figure it out as well.

Sourya

WARNING: Files paths have to be sorted.
WARNING: This was using a quite an old version of pomp -- some of the things will have to be adapted.

<params.R> model parameters
<init	.rda> I think these were initial conditions
<filter.study.setup.R> I setup the study
<fs.delta.chi.R> This file might have been how I put it on the cluster..
<msi_den_mod.R> <msi_den_mod.c> These are model files.
<casedata.csv> This is the simulated case data that I used to infer the results from.
<fs.delta.chi.rda> This is the result from the grid search
<prof.par1.R> This is how I generated likelihood profiles