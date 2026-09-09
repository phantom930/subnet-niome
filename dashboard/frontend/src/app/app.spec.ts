import { TestBed } from '@angular/core/testing';
import { provideRouter } from '@angular/router';
import { providePrimeNG } from 'primeng/config';

import { App } from './app';
import { AppPreset } from './theme/app-preset';

describe('App', () => {
  beforeEach(async () => {
    await TestBed.configureTestingModule({
      imports: [App],
      providers: [provideRouter([]), providePrimeNG({ theme: { preset: AppPreset } })],
    }).compileComponents();
  });

  it('should create the app', () => {
    const fixture = TestBed.createComponent(App);
    expect(fixture.componentInstance).toBeTruthy();
  });

  it('should render the toolbar shell and section navigation', async () => {
    const fixture = TestBed.createComponent(App);
    await fixture.whenStable();
    const compiled = fixture.nativeElement as HTMLElement;

    expect(compiled.querySelector('.p-toolbar')).toBeTruthy();
    expect(compiled.querySelectorAll('.p-button').length).toBeGreaterThan(0);

    const links = Array.from(compiled.querySelectorAll('.nav-link')).map((el) =>
      el.textContent?.trim(),
    );
    expect(links).toEqual(['Tasks', 'Mining']);
  });

  it('should toggle the dark mode class on the document root', async () => {
    const fixture = TestBed.createComponent(App);
    await fixture.whenStable();

    const root = document.documentElement;
    const before = root.classList.contains('app-dark');

    const toggle = fixture.nativeElement.querySelector('button') as HTMLButtonElement;
    toggle.click();
    await fixture.whenStable();

    expect(root.classList.contains('app-dark')).toBe(!before);
  });
});
