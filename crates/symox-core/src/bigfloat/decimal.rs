//! Certified decimal conversion with precision-sized scaled integers.
//!
//! Positive dyadic intervals are rounded outwards after each product/division.
//! Powers of ten use binary exponentiation; their exponent is stored separately
//! as an `Exponent` so a large decimal exponent does not allocate that many
//! digits.
use super::{BigFloat,Exponent};
use num_bigint::BigInt;
use num_integer::Integer;
use num_traits::{One,Signed,Zero};

#[derive(Clone,PartialEq,Eq)]
struct Positive {man:BigInt,exp:Exponent}
impl Positive {
    fn new(mut man:BigInt,exp:Exponent)->Self{
        if man.is_zero(){return Self{man,exp:Exponent::ZERO};}
        let zeros=man.trailing_zeros().unwrap();
        man>>=zeros as usize;
        Self{man,exp:exp+zeros}
    }
    fn round(&self,bits:u64,up:bool)->Self{
        let shift=self.man.bits().saturating_sub(bits.max(1));
        if shift==0{return self.clone();}
        let discarded=self.man.trailing_zeros().unwrap_or(shift)<shift;
        let mut man=&self.man>>shift as usize;
        if up&&discarded{man+=1u8;}
        Self::new(man,&self.exp+shift)
    }
    fn nearest(&self,bits:u64)->Self{
        let shift=self.man.bits().saturating_sub(bits.max(1));
        if shift==0{return self.clone();}
        let mut man=&self.man>>shift as usize;
        let half_bit=(&self.man>>(shift-1) as usize).is_odd();
        let above_half=self.man.trailing_zeros().unwrap_or(shift)<shift-1;
        if half_bit&&(above_half||man.is_odd()){man+=1u8;}
        Self::new(man,&self.exp+shift)
    }
    fn product(&self,other:&Self,bits:u64,up:bool)->Self{
        Self::new(&self.man*&other.man,&self.exp+&other.exp).round(bits,up)
    }
    fn quotient(&self,other:&Self,bits:u64,up:bool)->Self{
        let shift=(bits as i128+other.man.bits() as i128-self.man.bits() as i128+2).max(0) as usize;
        let (mut man,remainder)=(&self.man<<shift).div_rem(&other.man);
        if up&&!remainder.is_zero(){man+=1u8;}
        Self::new(man,&self.exp-&other.exp-shift).round(bits,up)
    }
    fn nearest_integer(&self)->BigInt{
        if self.exp>=0{return &self.man<<self.exp.to_usize().expect("integer exponent exceeds addressable storage");}
        let shift=match (-&self.exp).to_usize(){Some(shift) if shift as u64<=self.man.bits()=>shift,_=>return BigInt::zero()};
        let mut integer=&self.man>>shift;
        let half_bit=(&self.man>>(shift-1)).is_odd();
        // mpmath/SymPy decimal display resolves exact halves away from zero.
        // Binary construction instead uses nearest-even in nearest().
        if half_bit{integer+=1u8;}
        integer
    }
    fn bigfloat(self)->Result<BigFloat,String>{
        if self.exp.bits()>=u64::from(super::MAX_EXPONENT_BITS){
            return Err("decimal value exceeds the native binary exponent counter".into());
        }
        Ok(BigFloat{man:self.man,exp:self.exp})
    }
}

fn signed_round(man:BigInt,exp:Exponent,precision:u32)->Option<BigFloat>{
    let negative=man.is_negative();
    let result=Positive::new(man.abs(),exp).nearest(precision as u64).bigfloat().ok()?;
    Some(if negative{result.neg()}else{result})
}

pub(super) fn checked_product(a:&BigFloat,b:&BigFloat,precision:u32)->Option<BigFloat>{
    signed_round(&a.man*&b.man,&a.exp+&b.exp,precision)
}

pub(super) fn checked_sum(a:&BigFloat,b:&BigFloat,precision:u32)->Option<BigFloat>{
    if a.is_zero(){return checked_round(b,precision);}
    if b.is_zero(){return checked_round(a,precision);}
    let am=a.mag_exact();let bm=b.mag_exact();
    let (a,b,gap)=if am>=bm{(a,b,am-bm)}else{(b,a,bm-am)};
    if gap>precision as i128+4{
        let shift=(precision as i128+4-a.man.bits() as i128).max(0) as usize;
        let man=(&a.man<<(shift+2))+if b.is_negative(){-1i8}else{1i8};
        return signed_round(man,&a.exp-(shift+2),precision);
    }
    let exponent=a.exp.min_ref(&b.exp);
    let man=(&a.man<<(&a.exp-exponent).to_usize()?)
           +(&b.man<<(&b.exp-exponent).to_usize()?);
    signed_round(man,exponent.clone(),precision)
}

pub(super) fn checked_round(value:&BigFloat,precision:u32)->Option<BigFloat>{
    let result=Positive::new(value.man.abs(),value.exp.clone()).nearest(precision as u64).bigfloat().ok()?;
    Some(if value.is_negative(){result.neg()}else{result})
}

#[cfg(test)]
mod tests{
    use super::*;
    use num_rational::BigRational;

    #[test]
    fn large_decimal_powers_match_independent_exact_ratios(){
        let mantissa=BigInt::from(123456789u64);
        for precision in [1,2,3,24,53,167,500]{
            for exponent in [-5000i64,-1000,-130,0,130,1000,5000]{
                let power=num_traits::pow::Pow::pow(BigInt::from(10),exponent.unsigned_abs());
                for negative in [false,true]{
                    let n=if negative{-mantissa.clone()}else{mantissa.clone()};
                    let expected=if exponent>=0{BigFloat::from_int(&(n*&power)).round_to(precision)}
                                 else{BigFloat::from_ratio(&n,&power,precision)};
                    let literal=format!("{}123456789e{exponent}",if negative{"-"}else{""});
                    assert_eq!(from_decimal(&literal,precision).unwrap(),expected,"{literal} at {precision} bits");
                }
            }
        }
    }

    #[test]
    fn decimal_digits_have_independently_rounded_exact_values(){
        for mantissa in [1u64,3,5,15,125,123456789,987654321]{
            for exponent in [-2000i64,-200,-3,0,3,200,2000]{
                let value=BigFloat{man:mantissa.into(),exp:exponent.into()};
                for requested in [1u32,2,3,15,30,60]{
                    let (negative,text,decimal_exp)=digits(&value,requested);
                    assert!(!negative);assert_eq!(text.len(),requested as usize);
                    assert_ne!(text.as_bytes()[0],b'0');
                    let mut rational=value.to_ratio_exact();
                    let power=requested as i128-1-decimal_exp.small().unwrap();
                    let ten=num_traits::pow::Pow::pow(BigInt::from(10),power.unsigned_abs() as u64);
                    if power>=0{rational*=BigRational::from_integer(ten);}else{rational/=BigRational::from_integer(ten);}
                    let (mut integer,remainder)=rational.numer().div_rem(rational.denom());
                    let comparison=(&remainder<<1usize).cmp(rational.denom());
                    if comparison.is_ge(){integer+=1u8;}
                    assert_eq!(text,integer.to_string());
                }
            }
        }
    }

    #[test]
    fn integer_powers_match_independent_exact_ratios(){
        for mantissa in [-257i64,-31,-7,-3,-1,1,3,7,31,257]{
            for exponent in [-40i64,-1,0,1,40]{
                let value=BigFloat{man:mantissa.into(),exp:exponent.into()};
                for power in [-31i32,-3,-2,-1,0,1,2,3,31]{
                    let exact=value.to_ratio_exact().pow(power);
                    for precision in [1u32,3,53,167]{
                        let expected=BigFloat::from_ratio(exact.numer(),exact.denom(),precision);
                        assert_eq!(checked_integer_power(&value,power as i64,precision).unwrap(),expected);
                    }
                }
            }
        }
    }

    #[test]
    fn exponent_counter_overflow_is_held_instead_of_wrapping(){
        let value=BigFloat::from_decimal("1e999999999999999999",70).unwrap();
        let value=crate::Expr::num(crate::number::Number::mp_float(value,70));
        // Binary exponents are big integers of fewer than MAX_EXPONENT_BITS
        // bits. Each squaring doubles their size, so the bound is reached
        // after a few dozen steps.
        let mut square=crate::ops::mul2(&value,&value);
        while crate::ops::mul2(&square,&square).as_number().is_some(){
            square=crate::ops::mul2(&square,&square);
        }
        assert!(square.as_number().is_some());
        let cube=crate::ops::mul2(&square,&square);
        assert!(cube.is_mul());assert!(!cube.is_canonical());
        assert!(cube.has_wide_exponent());
        assert!(crate::evalf_mp::eval_mp(&cube,100).is_none());
        assert_eq!(crate::ops::add2(&square,&crate::ops::inv(&square)),square);
        let unit=crate::ops::div(&value,&value);
        assert_eq!(unit.as_number().unwrap().to_f64(),1.0);
    }

    #[test]
    fn trillion_digit_exponents_use_precision_sized_work(){
        for exponent in [-1_000_000_000_000i64,1_000_000_000_000]{
            let literal=format!("1e{exponent}");
            let value=from_decimal(&literal,53).unwrap();
            assert!(value.man.bits()<=53);
            let (_,text,decimal_exp)=digits(&value,15);
            assert_eq!(text,"100000000000000");assert_eq!(decimal_exp,i128::from(exponent));
        }
        // Decimal exponents beyond i64, as SymPy's Float('1e100000000000000000000').
        for exponent in [-100_000_000_000_000_000_000i128,100_000_000_000_000_000_000,123_456_789_012_345_678_901_234_567]{
            let value=from_decimal(&format!("1e{exponent}"),53).unwrap();
            let (_,text,decimal_exp)=digits(&value,15);
            assert_eq!(text,"100000000000000");assert_eq!(decimal_exp,exponent);
        }
        // Binary exponents beyond i128 are big integers, as mpmath's.
        for literal in ["1e10000000000000000000000000000000000000","1e1000000000000000000000000000000000000000000000000000","-1.5e-1000000000000000000000000000000000000000000000000000"]{
            let value=from_decimal(literal,53).unwrap();
            let (_,text,decimal_exp)=digits(&value,15);
            let expected:num_bigint::BigInt=literal.split('e').nth(1).unwrap().parse().unwrap();
            assert_eq!(text,if literal.starts_with('-'){"150000000000000"}else{"100000000000000"});
            assert_eq!(decimal_exp.to_bigint(),expected);
        }
        for exponent in [i64::MIN,i64::MAX]{
            let value=BigFloat{man:3.into(),exp:exponent.into()};
            assert_eq!(digits(&value,30).1.len(),30);
        }
    }
}

struct Interval {lo:Positive,hi:Positive}
impl Interval {
    fn point(man:BigInt,exp:Exponent,bits:u64)->Self{
        let point=Positive::new(man,exp);
        Self{lo:point.round(bits,false),hi:point.round(bits,true)}
    }
    fn product(&self,other:&Self,bits:u64)->Self{
        Self{lo:self.lo.product(&other.lo,bits,false),hi:self.hi.product(&other.hi,bits,true)}
    }
    fn quotient(&self,other:&Self,bits:u64)->Self{
        Self{lo:self.lo.quotient(&other.hi,bits,false),hi:self.hi.quotient(&other.lo,bits,true)}
    }
}

fn scale(man:&BigInt,exp:Exponent,decimal_power:&Exponent,bits:u64)->Interval{
    let mut power=Interval::point(BigInt::one(),Exponent::ZERO,bits);
    let mut base=Interval::point(BigInt::from(10),Exponent::ZERO,bits);
    let exponent=decimal_power.abs().to_bigint();
    let length=exponent.bits();
    for bit in 0..length{
        if exponent.bit(bit){power=power.product(&base,bits);}
        if bit+1<length{base=base.product(&base,bits);}
    }
    let value=Interval::point(man.clone(),exp,bits);
    if decimal_power.is_negative(){value.quotient(&power,bits)}else{value.product(&power,bits)}
}

fn initial_bits(target:u64,power:&Exponent)->u64{
    target+64+power.bits()
}

/// Integer powers have the same enclosure proof as decimal scaling. Positive
/// dyadic ties become exact once the working mantissa reaches p+1 bits;
/// reciprocals with an odd denominator cannot be exact binary rounding ties.
pub(super) fn checked_integer_power(value:&BigFloat,exponent:i64,precision:u32)->Option<BigFloat>{
    if exponent==0{return Some(BigFloat::one());}
    if value.is_zero(){return if exponent>0{Some(BigFloat::zero())}else{None};}
    let negative=value.is_negative()&&exponent.unsigned_abs()&1!=0;
    let mut bits=initial_bits(precision as u64,&Exponent::from(exponent));
    loop{
        let mut power=Interval::point(BigInt::one(),Exponent::ZERO,bits);
        let mut base=Interval::point(value.man.abs(),value.exp.clone(),bits);
        let mut count=exponent.unsigned_abs();
        while count!=0{
            if count&1!=0{power=power.product(&base,bits);}
            count>>=1;
            if count!=0{base=base.product(&base,bits);}
        }
        if exponent<0{power=Interval::point(BigInt::one(),Exponent::ZERO,bits).quotient(&power,bits);}
        let lo=power.lo.nearest(precision as u64);let hi=power.hi.nearest(precision as u64);
        if lo==hi{let value=lo.bigfloat().ok()?;return Some(if negative{value.neg()}else{value});}
        bits=bits.checked_mul(2)?;
    }
}

/// The exact path handles possible decimal/binary rounding ties with integers
/// whose sizes are bounded by the input mantissa and requested output size.
fn exact_scaled_integer(value:&BigFloat,power:i128)->BigInt{
    let five=num_traits::pow::Pow::pow(BigInt::from(5),u64::try_from(power.unsigned_abs()).expect("exact decimal scale exceeds addressable storage"));
    let exp=&value.exp+power;
    if power>=0{return Positive::new(value.man.abs()*five,exp).nearest_integer();}
    let mut numerator=value.man.abs();let mut denominator=five;
    if exp>=0{numerator<<=exp.to_usize().unwrap();}
    else{denominator<<=(-&exp).to_usize().unwrap();}
    let (mut quotient,remainder)=numerator.div_rem(&denominator);
    match (&remainder<<1usize).cmp(&denominator){
        std::cmp::Ordering::Greater=>quotient+=1u8,
        std::cmp::Ordering::Equal=>quotient+=1u8,
        _=>{},
    }
    quotient
}

pub(super) fn from_decimal(text:&str,precision:u32)->Result<BigFloat,String>{
    let text=text.trim();
    let (mantissa,exponent)=match text.find(['e','E']){
        Some(i)=>(&text[..i],text[i+1..].parse::<BigInt>().map_err(|_|"invalid decimal exponent")?),
        None=>(text,BigInt::zero()),
    };
    let (negative,mantissa)=match mantissa.as_bytes().first(){
        Some(b'-')=>(true,&mantissa[1..]),Some(b'+')=>(false,&mantissa[1..]),_=>(false,mantissa),
    };
    let mut digits=String::new();let mut fractional=None;
    for character in mantissa.chars(){
        if character=='.'&&fractional.is_none(){fractional=Some(0i128);}
        else if character.is_ascii_digit(){digits.push(character);if let Some(ref mut n)=fractional{*n+=1;}}
        else{return Err("invalid decimal literal".into());}
    }
    if digits.is_empty(){return Err("invalid decimal literal".into());}
    let power=Exponent::from(exponent)-fractional.unwrap_or(0);
    let mut mantissa:BigInt=digits.parse().map_err(|_|"invalid decimal literal")?;
    if mantissa.is_zero(){return Ok(BigFloat::zero());}
    if power.unsigned_abs()<=u128::from(mantissa.bits()+precision as u64+64){
        let power=power.small().unwrap();
        let ten=num_traits::pow::Pow::pow(BigInt::from(10),power.unsigned_abs() as u64);
        if negative{mantissa=-mantissa;}
        let result=if power>=0{BigFloat::from_int(&(mantissa*ten)).round_to(precision)}
                   else{BigFloat::from_ratio(&mantissa,&ten,precision)};
        return Ok(result);
    }
    let mut bits=initial_bits(precision as u64,&power);
    loop{
        let bounds=scale(&mantissa,Exponent::ZERO,&power,bits);
        let lo=bounds.lo.nearest(precision as u64);let hi=bounds.hi.nearest(precision as u64);
        if lo==hi{let value=lo.bigfloat()?;return Ok(if negative{value.neg()}else{value});}
        bits=bits.checked_mul(2).ok_or("decimal conversion working precision overflow")?;
    }
}

/// floor(n*log10(2)) up to one unit, without the 53-bit error of an f64
/// estimate at binary exponents beyond 2**53.
fn decimal_exponent_estimate(n:&Exponent)->Exponent{
    if n.unsigned_abs()<1<<50{return Exponent::from(((n.small().unwrap() as f64)*std::f64::consts::LOG10_2).floor() as i128);}
    // The product keeps 200 bits beyond the integer part.
    let wp=200+n.bits() as u32;
    let ratio=super::ln2(wp).div(&BigFloat::from_i64(10).ln(wp).unwrap(),wp).unwrap();
    Exponent::from(BigFloat::norm(n.to_bigint(),0).mul(&ratio,wp).to_int_trunc())
}

pub(super) fn digits(value:&BigFloat,digits:u32)->(bool,String,Exponent){
    let digits=digits.max(1);
    let magnitude=value.mag_exact();
    let mut exponent=decimal_exponent_estimate(&(&magnitude-1i64));
    loop{
        let power=Exponent::from(digits as i128-1)-&exponent;
        let exact=power.unsigned_abs()<=u128::from(value.man.bits()+4*digits as u64+64);
        let (integer,scaled_magnitude)=if exact{(exact_scaled_integer(value,power.small().unwrap()),None)}else{
            let mut bits=initial_bits(4*digits as u64,&power);
            loop{
                let bounds=scale(&value.man.abs(),value.exp.clone(),&power,bits);
                let lo=bounds.lo.nearest_integer();let hi=bounds.hi.nearest_integer();
                if lo==hi{break (lo,Some(&bounds.lo.exp+bounds.lo.man.bits()));}
                bits=bits.checked_mul(2).expect("decimal conversion working precision overflow");
            }
        };
        let text=integer.to_string();
        if text.len() as u64>digits as u64{exponent=exponent+(text.len() as i128-digits as i128);continue;}
        if (text.len() as u64)<digits as u64||integer.is_zero(){
            if integer.is_zero(){
                // A floating exponent estimate can be several hundred off
                // near i64 boundaries. Its only role is choosing a scale;
                // the final digit count and integer rounding are certified.
                let binary=scaled_magnitude.map(|m|m-1i64).unwrap_or_else(||&magnitude-1i64+(power.small().unwrap() as f64*std::f64::consts::LOG2_10).floor() as i128);
                let correction=decimal_exponent_estimate(&binary)-digits+1i64;
                exponent=exponent+correction.min(Exponent::from(-1i64));
            }else{exponent=exponent-(digits as i128-text.len() as i128);}
            continue;
        }
        return (value.is_negative(),text,exponent);
    }
}
