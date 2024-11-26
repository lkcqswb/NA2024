def to_base_beta_fraction(n, beta, precision=20):
    """将小数 n 转换为 beta 进制表示，保留指定的小数位数"""
    result = []
    for _ in range(precision):
        n *= beta
        digit = int(n)
        result.append(digit)
        n -= digit
    return result
def weighted_sum(fraction_digits, beta, p, result):

    sum_value = sum(fraction_digits[i] / (beta ** (i + 1)) for i in range(0,p))
    
    return sum_value - result
def get(fraction_digits,index):
    if(index>=len(fraction_digits)-1):
        return 0
    return fraction_digits[index]

def find_ab_for_beta(beta,p):
    for A in range(1, beta**p):
        for B in range(1, beta**p):
            while(A>=beta): A/=beta
            while(B>=beta): B/=beta
            if(A>=B): continue
                
            # 计算 A / B
            result = A / B
            
            # 转换成 beta 进制的小数部分
            fraction_digits = to_base_beta_fraction(result, beta,3*p)
            sum=0
            for i in range(2*p-1,len(fraction_digits)):
                 sum+=fraction_digits[i]/beta**(i-2*p+2)

            if(sum>=0.5):
                    fraction_digits[2*p-2]+=1
                    for i in range(2*p-1,len(fraction_digits)):
                         fraction_digits[i]=0
            sum=0
            for i in range(p,2*p-1):
                 sum+=fraction_digits[i]/beta**(i-p+1)
            if(sum>=0.5):
                    fraction_digits[p-1]+=1
                    
            if(abs(weighted_sum(fraction_digits,beta,p,result))/result>=0.5*beta**(1-p)):
                    print(f"A = {A}, B = {B}, A / B = {result:.5f}, Fraction in base-{beta}: {fraction_digits}")
                    print(str(p)+" yes "+str(beta))


# 使用函数

for p in range(1,10):
    for beta in range(50,100):
        print(beta)
        find_ab_for_beta(beta,p)
