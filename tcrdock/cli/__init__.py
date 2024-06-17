def normalise(name):
    return name.replace("_", "-")

context_settings = {"token_normalize_func": normalise}