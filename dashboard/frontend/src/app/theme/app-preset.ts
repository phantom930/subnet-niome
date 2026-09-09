import Aura from '@primeuix/themes/aura';
import { definePreset } from '@primeuix/themes';

/**
 * PrimeNG theme preset for the app.
 *
 * Built on Aura with the primary palette remapped to emerald. Edit the values
 * below to re-brand, or drop this file and pass `Aura` straight to
 * providePrimeNG() for the stock theme.
 */
export const AppPreset = definePreset(Aura, {
  semantic: {
    primary: {
      50: '{emerald.50}',
      100: '{emerald.100}',
      200: '{emerald.200}',
      300: '{emerald.300}',
      400: '{emerald.400}',
      500: '{emerald.500}',
      600: '{emerald.600}',
      700: '{emerald.700}',
      800: '{emerald.800}',
      900: '{emerald.900}',
      950: '{emerald.950}',
    },
  },
});
