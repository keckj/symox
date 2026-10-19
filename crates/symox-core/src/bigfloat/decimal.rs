//! Certified decimal conversion with precision-sized scaled integers.
//!
//! Positive dyadic intervals are rounded outwards after each product/division.
//! Powers of ten use binary exponentiation; their exponent is stored separately
//! in i128 so a large decimal exponent does not allocate that many digits.
use super::BigFloat;
use num_bigint::BigInt;
use num_integer::Integer;
use num_traits::{One,Signed,Zero};

#[derive(Clone,PartialEq,Eq)]
struct Positive {man:BigInt,exp:i128}
impl Positive {
    fn new(mut man:BigInt,mut exp:i128)->Self{
        if man.is_zero(){return Self{man,exp:0};}
        let zeros=man.trailing_zeros().unwrap();
        man>>=zeros as usize;exp+=zeros as i128;
        Self{man,exp}
    }
    fn round(&self,bits:u64,up:bool)->Self{
        let shift=self.man.bits().saturating_sub(bits.max(1));
        if shift==0{return self.clone();}
        let discarded=self.man.trailing_zeros().unwrap_or(shift)<shift;
        let mut man=&self.man>>shift as usize;
        if up&&discarded{man+=1u8;}
        Self::new(man,self.exp+shift as i128)
    }
    fn nearest(&self,bits:u64)->Self{
        let shift=self.man.bits().saturating_sub(bits.max(1));
        if shift==0{return self.clone();}
        let mut man=&self.man>>shift as usize;
        let half_bit=(&self.man>>(shift-1) as usize).is_odd();
        let above_half=self.man.trailing_zeros().unwrap_or(shift)<shift-1;
        if half_bit&&(above_half||man.is_odd()){man+=1u8;}
        Self::new(man,self.exp+shift as i128)
    }
    fn product(&self,other:&Self,bits:u64,up:bool)->Self{
        Self::new(&self.man*&other.man,self.exp+other.exp).round(bits,up)
    }
    fn quotient(&self,other:&Self,bits:u64,up:bool)->Self{
        let shift=(bits as i128+other.man.bits() as i128-self.man.bits() as i128+2).max(0) as usize;
        let (mut man,remainder)=(&self.man<<shift).div_rem(&other.man);
        if up&&!remainder.is_zero(){man+=1u8;}
        Self::new(man,self.exp-other.exp-shift as i128).round(bits,up)
    }
    fn nearest_integer(&self)->BigInt{
        if self.exp>=0{return &self.man<<usize::try_from(self.exp).expect("integer exponent exceeds addressable storage");}
        let shift=self.exp.unsigned_abs();
        if shift>self.man.bits() as u128{return BigInt::zero();}
        let shift=shift as usize;
        let mut integer=&self.man>>shift;
        let half_bit=(&self.man>>(shift-1)).is_odd();
        // mpmath/SymPy decimal display resolves exact halves away from zero.
        // Binary construction instead uses nearest-even in nearest().
        if half_bit{integer+=1u8;}
        integer
    }
    fn bigfloat(self)->Result<BigFloat,String>{
        let exp=i64::try_from(self.exp).map_err(|_|"decimal value exceeds the native binary exponent counter")?;
        Ok(BigFloat{man:self.man,exp})
    }
}

fn signed_round(man:BigInt,exp:i128,precision:u32)->Option<BigFloat>{
    let negative=man.is_negative();
    let result=Positive::new(man.abs(),exp).nearest(precision as u64).bigfloat().ok()?;
    Some(if negative{result.neg()}else{result})
}

pub(super) fn checked_product(a:&BigFloat,b:&BigFloat,precision:u32)->Option<BigFloat>{
    signed_round(&a.man*&b.man,a.exp as i128+b.exp as i128,precision)
}

pub(super) fn checked_sum(a:&BigFloat,b:&BigFloat,precision:u32)->Option<BigFloat>{
    if a.is_zero(){return checked_round(b,precision);}
    if b.is_zero(){return checked_round(a,precision);}
    let am=a.exp as i128+a.man.bits() as i128;let bm=b.exp as i128+b.man.bits() as i128;
    let (a,b,gap)=if am>=bm{(a,b,am-bm)}else{(b,a,bm-am)};
    if gap>precision as i128+4{
        let shift=(precision as i128+4-a.man.bits() as i128).max(0) as usize;
        let man=(&a.man<<(shift+2))+if b.is_negative(){-1i8}else{1i8};
        return signed_round(man,a.exp as i128-shift as i128-2,precision);
    }
    let exponent=a.exp.min(b.exp);
    let man=(&a.man<<usize::try_from(a.exp as i128-exponent as i128).ok()?)
           +(&b.man<<usize::try_from(b.exp as i128-exponent as i128).ok()?);
    signed_round(man,exponent as i128,precision)
}

pub(super) fn checked_quotient(a:&BigFloat,b:&BigFloat,precision:u32)->Option<BigFloat>{
    if b.is_zero(){return None;}
    if a.is_zero(){return Some(BigFloat::zero());}
    let shift=(precision as i128+3+b.man.bits() as i128-a.man.bits() as i128).max(0) as usize;
    let (quotient,remainder)=(&a.man.abs()<<shift).div_rem(&b.man.abs());
    let man=(quotient<<1usize)+u8::from(!remainder.is_zero());
    let man=if a.is_negative()!=b.is_negative(){-man}else{man};
    signed_round(man,a.exp as i128-b.exp as i128-shift as i128-1,precision)
}

pub(super) fn checked_round(value:&BigFloat,precision:u32)->Option<BigFloat>{
    let result=Positive::new(value.man.abs(),value.exp as i128).nearest(precision as u64).bigfloat().ok()?;
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
                let value=BigFloat{man:mantissa.into(),exp:exponent};
                for requested in [1u32,2,3,15,30,60]{
                    let (negative,text,decimal_exp)=digits(&value,requested);
                    assert!(!negative);assert_eq!(text.len(),requested as usize);
                    assert_ne!(text.as_bytes()[0],b'0');
                    let mut rational=value.to_ratio_exact();
                    let power=requested as i64-1-decimal_exp;
                    let ten=num_traits::pow::Pow::pow(BigInt::from(10),power.unsigned_abs());
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
                let value=BigFloat{man:mantissa.into(),exp:exponent};
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
        let square=crate::ops::mul2(&value,&value);
        assert!(square.as_number().is_some());
        let cube=crate::ops::mul2(&square,&value);
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
            assert_eq!(text,"100000000000000");assert_eq!(decimal_exp,exponent);
        }
        for exponent in [i64::MIN,i64::MAX]{
            let value=BigFloat{man:3.into(),exp:exponent};
            assert_eq!(digits(&value,30).1.len(),30);
        }
    }
}

struct Interval {lo:Positive,hi:Positive}
impl Interval {
    fn point(man:BigInt,exp:i128,bits:u64)->Self{
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

fn scale(man:&BigInt,exp:i128,decimal_power:i64,bits:u64)->Interval{
    let mut power=Interval::point(BigInt::one(),0,bits);
    let mut base=Interval::point(BigInt::from(10),0,bits);
    let mut exponent=decimal_power.unsigned_abs();
    while exponent!=0{
        if exponent&1!=0{power=power.product(&base,bits);}
        exponent>>=1;
        if exponent!=0{base=base.product(&base,bits);}
    }
    let value=Interval::point(man.clone(),exp,bits);
    if decimal_power<0{value.quotient(&power,bits)}else{value.product(&power,bits)}
}

fn initial_bits(target:u64,power:i64)->u64{
    target+64+(64-power.unsigned_abs().leading_zeros()) as u64
}

/// Integer powers have the same enclosure proof as decimal scaling. Positive
/// dyadic ties become exact once the working mantissa reaches p+1 bits;
/// reciprocals with an odd denominator cannot be exact binary rounding ties.
pub(super) fn checked_integer_power(value:&BigFloat,exponent:i64,precision:u32)->Option<BigFloat>{
    if exponent==0{return Some(BigFloat::one());}
    if value.is_zero(){return if exponent>0{Some(BigFloat::zero())}else{None};}
    let negative=value.is_negative()&&exponent.unsigned_abs()&1!=0;
    let mut bits=initial_bits(precision as u64,exponent);
    loop{
        let mut power=Interval::point(BigInt::one(),0,bits);
        let mut base=Interval::point(value.man.abs(),value.exp as i128,bits);
        let mut count=exponent.unsigned_abs();
        while count!=0{
            if count&1!=0{power=power.product(&base,bits);}
            count>>=1;
            if count!=0{base=base.product(&base,bits);}
        }
        if exponent<0{power=Interval::point(BigInt::one(),0,bits).quotient(&power,bits);}
        let lo=power.lo.nearest(precision as u64);let hi=power.hi.nearest(precision as u64);
        if lo==hi{let value=lo.bigfloat().ok()?;return Some(if negative{value.neg()}else{value});}
        bits=bits.checked_mul(2)?;
    }
}

/// The exact path handles possible decimal/binary rounding ties with integers
/// whose sizes are bounded by the input mantissa and requested output size.
fn exact_scaled_integer(value:&BigFloat,power:i64)->BigInt{
    let five=num_traits::pow::Pow::pow(BigInt::from(5),power.unsigned_abs());
    let exp=value.exp as i128+power as i128;
    if power>=0{return Positive::new(value.man.abs()*five,exp).nearest_integer();}
    let mut numerator=value.man.abs();let mut denominator=five;
    if exp>=0{numerator<<=usize::try_from(exp).unwrap();}
    else{denominator<<=usize::try_from(exp.unsigned_abs()).unwrap();}
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
        Some(i)=>(&text[..i],text[i+1..].parse::<i64>().map_err(|_|"invalid decimal exponent")?),
        None=>(text,0),
    };
    let (negative,mantissa)=match mantissa.as_bytes().first(){
        Some(b'-')=>(true,&mantissa[1..]),Some(b'+')=>(false,&mantissa[1..]),_=>(false,mantissa),
    };
    let mut digits=String::new();let mut fractional=None;
    for character in mantissa.chars(){
        if character=='.'&&fractional.is_none(){fractional=Some(0i64);}
        else if character.is_ascii_digit(){digits.push(character);if let Some(ref mut n)=fractional{*n+=1;}}
        else{return Err("invalid decimal literal".into());}
    }
    if digits.is_empty(){return Err("invalid decimal literal".into());}
    let power=exponent.checked_sub(fractional.unwrap_or(0)).ok_or("decimal exponent overflow")?;
    let mut mantissa:BigInt=digits.parse().map_err(|_|"invalid decimal literal")?;
    if mantissa.is_zero(){return Ok(BigFloat::zero());}
    if power.unsigned_abs()<=mantissa.bits()+precision as u64+64{
        let ten=num_traits::pow::Pow::pow(BigInt::from(10),power.unsigned_abs());
        if negative{mantissa=-mantissa;}
        let result=if power>=0{BigFloat::from_int(&(mantissa*ten)).round_to(precision)}
                   else{BigFloat::from_ratio(&mantissa,&ten,precision)};
        return Ok(result);
    }
    let mut bits=initial_bits(precision as u64,power);
    loop{
        let bounds=scale(&mantissa,0,power,bits);
        let lo=bounds.lo.nearest(precision as u64);let hi=bounds.hi.nearest(precision as u64);
        if lo==hi{let value=lo.bigfloat()?;return Ok(if negative{value.neg()}else{value});}
        bits=bits.checked_mul(2).ok_or("decimal conversion working precision overflow")?;
    }
}

pub(super) fn digits(value:&BigFloat,digits:u32)->(bool,String,i64){
    let digits=digits.max(1);
    let magnitude=value.exp as i128+value.man.bits() as i128;
    let mut exponent=(((magnitude-1) as f64)*std::f64::consts::LOG10_2).floor() as i64;
    loop{
        let power=digits as i64-1-exponent;
        let exact=power.unsigned_abs()<=value.man.bits()+4*digits as u64+64;
        let (integer,scaled_magnitude)=if exact{(exact_scaled_integer(value,power),None)}else{
            let mut bits=initial_bits(4*digits as u64,power);
            loop{
                let bounds=scale(&value.man.abs(),value.exp as i128,power,bits);
                let lo=bounds.lo.nearest_integer();let hi=bounds.hi.nearest_integer();
                if lo==hi{break (lo,Some(bounds.lo.exp+bounds.lo.man.bits() as i128));}
                bits=bits.checked_mul(2).expect("decimal conversion working precision overflow");
            }
        };
        let text=integer.to_string();
        if text.len() as u64>digits as u64{exponent+=text.len() as i64-digits as i64;continue;}
        if (text.len() as u64)<digits as u64||integer.is_zero(){
            if integer.is_zero(){
                // A floating exponent estimate can be several hundred off
                // near i64 boundaries. Its only role is choosing a scale;
                // the final digit count and integer rounding are certified.
                let binary=scaled_magnitude.map(|m|m-1).unwrap_or_else(||magnitude-1+(power as f64*std::f64::consts::LOG2_10).floor() as i128);
                let correction=((binary as f64)*std::f64::consts::LOG10_2).floor() as i64-digits as i64+1;
                exponent+=correction.min(-1);
            }else{exponent-=digits as i64-text.len() as i64;}
            continue;
        }
        return (value.is_negative(),text,exponent);
    }
}
