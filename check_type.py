def check_type(val):
    try:
      ans=int(val)
      return ans
    except ValueError:
      try:
        ans=float(val)
        return ans
      except ValueError:
        if val.lower()=="false":
          return False
        elif val.lower()=="true":
          return True
        else:
          return val


