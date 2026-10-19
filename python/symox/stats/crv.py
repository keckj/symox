"""Continuous domains, distributions and probability spaces."""
from ._domains import ContinuousDomain,SingleContinuousDomain,ProductContinuousDomain,ConditionalContinuousDomain
from ._spaces import ContinuousDistribution,SingleContinuousDistribution,ContinuousPSpace,SingleContinuousPSpace

def reduce_rational_inequalities_wrap(condition,var):
    from symox import solve_univariate_inequality
    return solve_univariate_inequality(condition,var,relational=False)

__all__=['ContinuousDomain','SingleContinuousDomain','ProductContinuousDomain','ConditionalContinuousDomain',
 'ContinuousDistribution','SingleContinuousDistribution','ContinuousPSpace','SingleContinuousPSpace','reduce_rational_inequalities_wrap']
