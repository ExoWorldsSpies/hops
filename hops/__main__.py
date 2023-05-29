
if __name__ == '__main__':
    import hops
    import traceback
    try:
        hops.run_app()
    except:
        print('\n************\n\n')
        traceback.print_exc()
        x = input('Press enter to exit.\n')
        raise