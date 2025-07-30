def rotl32(x, n):
    return ((x << n) & 0xffffffff) | (x >> (32 - n))

def quarter_round(a, b, c, d):
    a = (a + b) & 0xffffffff
    d ^= a
    d = rotl32(d, 16)
    c = (c + d) & 0xffffffff
    b ^= c
    b = rotl32(b, 12)
    a = (a + b) & 0xffffffff
    d ^= a
    d = rotl32(d, 8)
    c = (c + d) & 0xffffffff
    b ^= c
    b = rotl32(b, 7)
    return a, b, c, d

def chacha20_block(key, counter, nonce):
    # constant = [0, 0, 0, 0]
    state = [0, 0, 0, 0] + key + [counter] + nonce
    print("Initial state:", [hex(x) for x in state])
    working_state = state[:]
    # working_state[0], working_state[4], working_state[8], working_state[12] = quarter_round(
    #         working_state[0], working_state[4], working_state[8], working_state[12])
    # working_state[1], working_state[5], working_state[9], working_state[13] = quarter_round(
    #         working_state[1], working_state[5], working_state[9], working_state[13])
    # working_state[2], working_state[6], working_state[10], working_state[14] = quarter_round(
    #     working_state[2], working_state[6], working_state[10], working_state[14])
    # working_state[3], working_state[7], working_state[11], working_state[15] = quarter_round(
    #     working_state[3], working_state[7], working_state[11], working_state[15])
    for _ in range(10):
        # Column rounds
        working_state[0], working_state[4], working_state[8], working_state[12] = quarter_round(
            working_state[0], working_state[4], working_state[8], working_state[12])
        working_state[1], working_state[5], working_state[9], working_state[13] = quarter_round(
            working_state[1], working_state[5], working_state[9], working_state[13])
        working_state[2], working_state[6], working_state[10], working_state[14] = quarter_round(
            working_state[2], working_state[6], working_state[10], working_state[14])
        working_state[3], working_state[7], working_state[11], working_state[15] = quarter_round(
            working_state[3], working_state[7], working_state[11], working_state[15])
        # Diagonal rounds
        working_state[0], working_state[5], working_state[10], working_state[15] = quarter_round(
            working_state[0], working_state[5], working_state[10], working_state[15])
        working_state[1], working_state[6], working_state[11], working_state[12] = quarter_round(
            working_state[1], working_state[6], working_state[11], working_state[12])
        working_state[2], working_state[7], working_state[8], working_state[13] = quarter_round(
            working_state[2], working_state[7], working_state[8], working_state[13])
        working_state[3], working_state[4], working_state[9], working_state[14] = quarter_round(
            working_state[3], working_state[4], working_state[9], working_state[14])
    output = [(working_state[i] + state[i]) & 0xffffffff for i in range(16)]
    # return working_state
    return output

# Ví dụ sử dụng
if __name__ == "__main__":
    key = [0] * 8  # 256-bit key (8 x 32-bit)
    counter = 1
    nonce = [0, 0, 0]  # 96-bit nonce (3 x 32-bit)
    block = chacha20_block(key, counter, nonce)
    print("Key:", [hex(x) for x in key])
    print("Counter:", counter)
    print("Nonce:", [hex(x) for x in nonce])
    print("ChaCha20 block với constant = 0:")
    for x in block:
        print(hex(x), end=" ")