* Computes Tr(gamma(mu1) ... gamma(mu14)) in four dimensions.
* Run with: form trace_14_gammas.frm

#-

Indices mu1,...,mu17;

Local trace14 =
    g_(1,mu1,mu2,mu3,mu4,mu5,mu6,mu7,mu8,mu9,mu10,mu11,mu12,mu13,mu14);

trace4,1;
.sort

Print +S;
.end
