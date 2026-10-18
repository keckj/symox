//! Convergent series with explicit domain and iteration guards.
use super::{C64,cr,gamma};

fn nonpositive_integer(z:C64)->Option<usize>{
    if z.im==0.0 && z.re<=0.0 && z.re==z.re.floor() && z.re>=-2048.0 {Some((-z.re) as usize)}else{None}
}

pub(super) fn appellf1(a:&[C64])->Option<C64>{
    let (aa,b1,b2,c,x,y)=(a[0],a[1],a[2],a[3],a[4],a[5]);
    let degree=nonpositive_integer(aa).or_else(||Some(nonpositive_integer(b1)?+nonpositive_integer(b2)?));
    if degree.is_none() && (x.norm()>=0.98 || y.norm()>=0.98){return None;}
    let limit=degree.unwrap_or(2048);
    if limit>2048{return None;}
    let mut u=vec![cr(1.0)];let mut v=vec![cr(1.0)];let mut ratio=cr(1.0);
    let mut result=cr(1.0);let mut small=0;
    for k in 1..=limit {
        let prev=k as f64-1.0;let denominator=c+prev;
        if denominator.norm()==0.0{return None;}
        ratio*= (aa+prev)/denominator;
        u.push(u[k-1]*(b1+prev)*x/(k as f64));v.push(v[k-1]*(b2+prev)*y/(k as f64));
        let term=ratio*(0..=k).map(|j|u[j]*v[k-j]).sum::<C64>();
        result+=term;
        if !result.re.is_finite() || !result.im.is_finite(){return None;}
        if term.norm()<2e-16*result.norm().max(1.0){small+=1;}else{small=0;}
        if degree.is_none() && small>=8 && k>16{return Some(result);}
    }
    if degree.is_some(){Some(result)}else{None}
}

fn gamma_q(a:f64,x:f64)->Option<f64>{
    if a==0.0{return Some(0.0);}
    if x==0.0{return Some(1.0);}
    let prefactor=(a*x.ln()-x-gamma::loggamma_c(cr(a))?.re).exp();
    if x<a+1.0{
        let mut term=1.0/a;let mut sum=term;
        for k in 1..10000{
            term*=x/(a+k as f64);sum+=term;
            if term.abs()<2e-16*sum.abs(){return Some((1.0-prefactor*sum).clamp(0.0,1.0));}
        }
    }else{
        let tiny=1e-280;let mut b=x+1.0-a;let mut c=1.0/tiny;let mut d=1.0/b;let mut h=d;
        for k in 1..10000{
            let an=-(k as f64)*(k as f64-a);b+=2.0;
            d=an*d+b;if d.abs()<tiny{d=tiny;}
            c=b+an/c;if c.abs()<tiny{c=tiny;}
            d=1.0/d;let delta=d*c;h*=delta;
            if (delta-1.0).abs()<2e-15{return Some((prefactor*h).clamp(0.0,1.0));}
        }
    }
    None
}

pub(super) fn marcumq(a:&[C64])->Option<C64>{
    if a.iter().any(|z|z.im!=0.0 || !z.re.is_finite() || z.re<0.0){return None;}
    let (m,x,y)=(a[0].re,a[1].re,a[2].re);
    let lambda=x*x/2.0;let t=y*y/2.0;
    if lambda>4096.0{return None;}
    if lambda==0.0{return Some(cr(if m==0.0{0.0}else{gamma_q(m,t)?}));}
    let end=(lambda+14.0*(lambda+1.0).sqrt()+32.0).ceil() as usize;
    let mut sum=0.0;
    for k in 0..=end{
        let log_weight=-lambda+(k as f64)*lambda.ln()-gamma::loggamma_c(cr(k as f64+1.0))?.re;
        let weight=log_weight.exp();
        if weight>0.0{sum+=weight*gamma_q(m+k as f64,t)?;}
    }
    Some(cr(sum.clamp(0.0,1.0)))
}
